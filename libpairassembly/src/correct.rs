use std::{fmt::Display, sync::LazyLock};

use crate::{
    PairOverlap, Result,
    assembler::HasPairOverlap,
    errors::{ConversionError, CorrectionError::ConsensusLengthMismatch},
    merge::MergedConsensus,
    overlap::{HasOrientedPairSlices, OverlapBounds, private::Sealed},
    prelude::utils::encode_fastq_quality_scores_in_place,
    read::{OwnedReadPair, OwnedSequenceRead},
};

const MIN_EFFECTIVE_PHRED_INPUT: u8 = 0;
const MAX_EFFECTIVE_PHRED_INPUT: u8 = 41;
const MAX_CORRECTED_PHRED_OUTPUT: u8 = 40;
const QUALITY_TABLE_LEN: usize =
    (MAX_EFFECTIVE_PHRED_INPUT - MIN_EFFECTIVE_PHRED_INPUT + 1) as usize;

// Shared lookup tables for the correction kernel. These precompute error probabilities and the
// posterior-style corrected quality outputs for matching and mismatching overlap columns so the hot
// loop does not have to redo transcendental math at every base.
//
// As of the current Rust toolchain, the floating-point operations needed to build these tables
// (`powf`, `log10`) are not available in const evaluation, so these remain one-time initialized via
// `LazyLock` for now.
static CORRECTION_TABLES: LazyLock<CorrectionTables> = LazyLock::new(CorrectionTables::build);

/// Configuration for correction behavior.
///
/// This type is non-exhaustive so additional correction policies can be added later without
/// requiring a breaking change.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CorrectionParams {
    /// Maximum corrected quality score the kernel is allowed to emit.
    pub max_output_qual: u8,

    /// If true, correction updates overlap qualities only and leaves called bases unchanged.
    pub quality_only: bool,
}

impl Default for CorrectionParams {
    fn default() -> Self {
        Self {
            max_output_qual: MAX_CORRECTED_PHRED_OUTPUT,
            quality_only: false,
        }
    }
}

impl CorrectionParams {
    #[must_use]
    pub fn with_max_output_qual(self, max_output_qual: u8) -> Self {
        Self {
            max_output_qual,
            ..self
        }
    }

    #[must_use]
    pub fn quality_only(self) -> Self {
        Self {
            quality_only: true,
            ..self
        }
    }
}

/// Aligned overlap-local input window shared by correction kernels.
#[derive(Debug, Clone)]
struct CorrectionWindow<'a> {
    fwd_seq: &'a [u8],
    fwd_qual: &'a [u8],
    rev_seq: &'a [u8],
    rev_qual: &'a [u8],
}

/// Corrected consensus record emitted after applying overlap-based quality correction.
///
/// Internally, quality bytes are numeric scores. User-facing record egress encodes them as FASTQ
/// ASCII quality bytes while consuming the owned buffer.
#[derive(Debug, Clone)]
pub struct CorrectedMergedRead {
    id: String,
    seq: Vec<u8>,
    quality_scores: Vec<u8>,
}

/// Corrected oriented pair retained by staged assembler contexts.
///
/// This remains in score-space overlap orientation so corrected validation and merge stages can
/// reuse the corrected buffers directly.
#[derive(Debug, Clone)]
pub(crate) struct CorrectedOrientedPair {
    id: String,
    fwd_seq: Vec<u8>,
    fwd_quality_scores: Vec<u8>,
    rev_seq_rc: Vec<u8>,
    rev_quality_scores_rc: Vec<u8>,
    overlap_bounds: OverlapBounds,
}

/// Applies overlap-based correction to pair and merged-read slices.
pub(crate) struct OverlapCorrector {
    params: CorrectionParams,
    tables: &'static CorrectionTables,
}

impl<'a> CorrectionWindow<'a> {
    #[must_use]
    fn new(fwd_seq: &'a [u8], fwd_qual: &'a [u8], rev_seq: &'a [u8], rev_qual: &'a [u8]) -> Self {
        debug_assert_eq!(fwd_seq.len(), fwd_qual.len());
        debug_assert_eq!(rev_seq.len(), rev_qual.len());
        debug_assert_eq!(fwd_seq.len(), rev_seq.len());

        Self {
            fwd_seq,
            fwd_qual,
            rev_seq,
            rev_qual,
        }
    }

    #[must_use]
    fn from_overlap<'pair>(overlap: &'a PairOverlap<'pair>) -> Self {
        let (fwd_seq, rev_seq) = overlap.overlap_windows();
        let (fwd_qual, rev_qual) = overlap.overlap_quality_windows();

        Self::new(fwd_seq, fwd_qual, rev_seq, rev_qual)
    }

    #[must_use]
    fn len(&self) -> usize {
        self.fwd_seq.len()
    }

    #[must_use]
    fn forward_sequence(&self) -> &[u8] {
        self.fwd_seq
    }

    #[must_use]
    fn forward_qualities(&self) -> &[u8] {
        self.fwd_qual
    }

    #[must_use]
    fn reverse_sequence(&self) -> &[u8] {
        self.rev_seq
    }

    #[must_use]
    fn reverse_qualities(&self) -> &[u8] {
        self.rev_qual
    }
}

impl OverlapCorrector {
    pub(crate) fn new(params: CorrectionParams) -> Self {
        Self {
            params,
            tables: &CORRECTION_TABLES,
        }
    }

    pub(crate) fn correct_pair_overlap<T>(&self, target: &T) -> Result<CorrectedOrientedPair>
    where
        T: HasPairOverlap + ?Sized,
    {
        target.validate_overlap_bounds()?;

        let slices = target.pair_slices()?;
        let overlap_bounds = target.overlap_bounds()?;

        let mut fwd_seq = slices.forward_sequence().to_vec();
        let mut fwd_quality_scores = slices.forward_quality_score_bytes().to_vec();
        let mut rev_seq_rc = slices.reverse_sequence_rc().to_vec();
        let mut rev_quality_scores_rc = slices.reverse_quality_score_bytes_rc().to_vec();

        let fwd_range = overlap_bounds.forward_range();
        let rev_range = overlap_bounds.reverse_range();
        let (fwd_window_seq, rev_window_seq) = target.overlap_windows()?;
        let (fwd_window_qual, rev_window_qual) = target.overlap_quality_windows()?;
        let window = CorrectionWindow::new(
            fwd_window_seq,
            fwd_window_qual,
            rev_window_seq,
            rev_window_qual,
        );

        self.correct_oriented_pair_overlap(
            &window,
            &mut fwd_seq[fwd_range.clone()],
            &mut fwd_quality_scores[fwd_range],
            &mut rev_seq_rc[rev_range.clone()],
            &mut rev_quality_scores_rc[rev_range],
        );

        Ok(CorrectedOrientedPair {
            id: slices.pair_id().to_string(),
            fwd_seq,
            fwd_quality_scores,
            rev_seq_rc,
            rev_quality_scores_rc,
            overlap_bounds,
        })
    }

    pub(crate) fn correct_merged_consensus(
        &self,
        consensus: MergedConsensus,
        overlap: &PairOverlap<'_>,
    ) -> Result<CorrectedMergedRead> {
        let window = CorrectionWindow::from_overlap(overlap);
        let overlap_start = consensus.left_overhang_len();
        let MergedConsensus {
            id,
            sequence,
            quality_scores,
            ..
        } = consensus;

        self.correct_merged_parts(id, sequence, quality_scores, overlap_start, &window)
    }

    fn correct_merged_parts(
        &self,
        id: String,
        mut sequence: Vec<u8>,
        mut quality_scores: Vec<u8>,
        overlap_start: usize,
        window: &CorrectionWindow<'_>,
    ) -> Result<CorrectedMergedRead> {
        let overlap_end = overlap_start + window.len();

        self.correct_merged_overlap(
            window,
            &mut sequence[overlap_start..overlap_end],
            &mut quality_scores[overlap_start..overlap_end],
        );

        CorrectedMergedRead::try_new(id, sequence, quality_scores)
    }

    fn correct_merged_overlap(
        &self,
        window: &CorrectionWindow<'_>,
        seq_overlap: &mut [u8],
        qual_overlap: &mut [u8],
    ) {
        debug_assert_eq!(seq_overlap.len(), window.len());
        debug_assert_eq!(qual_overlap.len(), window.len());

        for idx in 0..window.len() {
            let (corrected_base, corrected_qual) = self.correct_column(window, idx);

            if !self.params.quality_only {
                seq_overlap[idx] = corrected_base;
            }
            qual_overlap[idx] = corrected_qual;
        }
    }

    fn correct_oriented_pair_overlap(
        &self,
        window: &CorrectionWindow<'_>,
        fwd_seq_overlap: &mut [u8],
        fwd_qual_overlap: &mut [u8],
        rev_seq_overlap_rc: &mut [u8],
        rev_qual_overlap_rc: &mut [u8],
    ) {
        debug_assert_eq!(fwd_seq_overlap.len(), window.len());
        debug_assert_eq!(fwd_qual_overlap.len(), window.len());
        debug_assert_eq!(rev_seq_overlap_rc.len(), window.len());
        debug_assert_eq!(rev_qual_overlap_rc.len(), window.len());

        for idx in 0..window.len() {
            let (corrected_base, corrected_qual) = self.correct_column(window, idx);

            let write_base = if self.params.quality_only {
                window.forward_sequence()[idx]
            } else {
                corrected_base
            };

            fwd_seq_overlap[idx] = write_base;
            fwd_qual_overlap[idx] = corrected_qual;
            rev_seq_overlap_rc[idx] = write_base;
            rev_qual_overlap_rc[idx] = corrected_qual;
        }
    }

    fn correct_column(&self, window: &CorrectionWindow<'_>, idx: usize) -> (u8, u8) {
        self.tables.correct_overlap_column(
            window.forward_sequence()[idx],
            window.reverse_sequence()[idx],
            window.forward_qualities()[idx],
            window.reverse_qualities()[idx],
            self.params,
        )
    }
}

impl CorrectedOrientedPair {
    #[inline]
    pub(crate) fn overlap_bounds(&self) -> OverlapBounds {
        self.overlap_bounds
    }

    /// Consume corrected pair slices into merged consensus buffers.
    ///
    /// The corrected forward sequence and quality buffers already contain the merged left overhang
    /// and corrected overlap in the right orientation. Merging a corrected pair can therefore reuse
    /// those buffers by truncating any forward-only suffix and appending the reverse mate's right
    /// overhang from the reverse-complemented slices.
    ///
    /// # Errors
    ///
    /// Returns an error if the retained overlap bounds are inconsistent with the corrected slices
    /// buffers, or if the resulting consensus violates sequence/quality length invariants.
    pub(crate) fn to_merged_consensus(&self) -> Result<MergedConsensus> {
        self.validate_overlap_bounds(self.overlap_bounds)?;

        let fwd_start = self.overlap_bounds.fwd_start_offset();
        let fwd_end = self.overlap_bounds.forward_range().end;
        let rev_end = self.overlap_bounds.reverse_range().end;

        let mut sequence = Vec::with_capacity(fwd_end + self.rev_seq_rc.len() - rev_end);
        sequence.extend_from_slice(&self.fwd_seq[..fwd_end]);
        sequence.extend_from_slice(&self.rev_seq_rc[rev_end..]);

        let mut quality_scores =
            Vec::with_capacity(fwd_end + self.rev_quality_scores_rc.len() - rev_end);
        quality_scores.extend_from_slice(&self.fwd_quality_scores[..fwd_end]);
        quality_scores.extend_from_slice(&self.rev_quality_scores_rc[rev_end..]);

        MergedConsensus::try_new(self.id.clone(), sequence, quality_scores, fwd_start)
    }
}

impl Sealed for CorrectedOrientedPair {}

impl HasOrientedPairSlices for CorrectedOrientedPair {
    fn pair_id(&self) -> &str {
        &self.id
    }

    fn forward_sequence(&self) -> &[u8] {
        &self.fwd_seq
    }

    fn forward_quality_score_bytes(&self) -> &[u8] {
        &self.fwd_quality_scores
    }

    fn reverse_sequence_rc(&self) -> &[u8] {
        &self.rev_seq_rc
    }

    fn reverse_quality_score_bytes_rc(&self) -> &[u8] {
        &self.rev_quality_scores_rc
    }
}

impl CorrectedMergedRead {
    /// Construct a corrected merged read from checked record parts.
    ///
    /// # Errors
    ///
    /// Returns an error if sequence and quality lengths differ.
    pub(crate) fn try_new(id: String, seq: Vec<u8>, quality_scores: Vec<u8>) -> Result<Self> {
        if seq.len() != quality_scores.len() {
            return Err(ConsensusLengthMismatch {
                seq_len: seq.len(),
                qual_len: quality_scores.len(),
            }
            .into());
        }

        Ok(Self {
            id,
            seq,
            quality_scores,
        })
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn sequence_bytes(&self) -> &[u8] {
        self.seq.as_slice()
    }

    #[must_use]
    pub fn sequence_owned(self) -> Vec<u8> {
        self.seq
    }

    #[must_use]
    pub fn quality_score_bytes(&self) -> &[u8] {
        self.quality_scores.as_slice()
    }

    #[must_use]
    pub fn quality_score_bytes_owned(self) -> Vec<u8> {
        self.quality_scores
    }

    /// Return ASCII-encoded FASTQ quality bytes for the corrected merged consensus.
    #[must_use]
    pub fn to_quality_ascii_bytes(&self) -> Vec<u8> {
        let mut quality_ascii = self.quality_scores.clone();
        encode_fastq_quality_scores_in_place(&mut quality_ascii);
        quality_ascii
    }

    /// Convert corrected merged output into a user record value.
    ///
    /// This is a boundary conversion API. Identity-shaped `into_*` methods are
    /// intentionally omitted to keep `into_*` naming reserved for meaningful
    /// representation changes.
    ///
    /// # Errors
    ///
    /// Returns an error if the target record type cannot be constructed
    /// from the corrected parts.
    pub fn into_record<T>(self) -> Result<T>
    where
        T: TryFrom<OwnedSequenceRead>,
        T::Error: Display,
    {
        let read = OwnedSequenceRead::try_from(self)?;
        T::try_from(read).map_err(|err| ConversionError::RecordConstruction(err.to_string()).into())
    }
}

impl TryFrom<CorrectedOrientedPair> for OwnedReadPair {
    type Error = crate::Error;

    fn try_from(mut corrected_pair: CorrectedOrientedPair) -> Result<Self> {
        let mut fwd_quality_ascii = corrected_pair.fwd_quality_scores;
        encode_fastq_quality_scores_in_place(&mut fwd_quality_ascii);

        corrected_pair.rev_seq_rc.reverse();
        for base in &mut corrected_pair.rev_seq_rc {
            *base = complement_base(*base);
        }

        let mut rev_quality_ascii = corrected_pair.rev_quality_scores_rc;
        rev_quality_ascii.reverse();
        encode_fastq_quality_scores_in_place(&mut rev_quality_ascii);

        OwnedReadPair::builder()
            .id(corrected_pair.id)
            .forward(corrected_pair.fwd_seq, fwd_quality_ascii)
            .reverse(corrected_pair.rev_seq_rc, rev_quality_ascii)
            .build()
    }
}

impl TryFrom<CorrectedMergedRead> for OwnedSequenceRead {
    type Error = crate::Error;

    fn try_from(read: CorrectedMergedRead) -> Result<Self> {
        let mut quality_ascii = read.quality_scores;
        encode_fastq_quality_scores_in_place(&mut quality_ascii);

        Self::try_from_ascii_bytes(read.id, read.seq, quality_ascii)
    }
}

impl TryFrom<MergedConsensus> for CorrectedMergedRead {
    type Error = crate::Error;

    fn try_from(consensus: MergedConsensus) -> Result<Self> {
        let MergedConsensus {
            id,
            sequence,
            quality_scores,
            ..
        } = consensus;

        Self::try_new(id, sequence, quality_scores)
    }
}

#[inline]
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

#[inline]
fn mismatch_error_probability(fwd_error: f64, rev_error: f64) -> f64 {
    ((fwd_error * rev_error) / 3.0)
        / ((1.0 - fwd_error) * (1.0 - rev_error) + 4.0 * (fwd_error * rev_error) / 3.0)
}

#[inline]
fn match_error_probability(fwd_error: f64, rev_error: f64) -> f64 {
    (fwd_error * (1.0 - rev_error / 3.0))
        / (fwd_error + rev_error - 4.0 * (fwd_error * rev_error) / 3.0)
}

#[derive(Debug)]
struct CorrectionTables {
    match_qual: [[u8; QUALITY_TABLE_LEN]; QUALITY_TABLE_LEN],
    mismatch_qual: [[u8; QUALITY_TABLE_LEN]; QUALITY_TABLE_LEN],
}

impl CorrectionTables {
    fn build() -> Self {
        let mut error_prob = [0.0f64; QUALITY_TABLE_LEN];
        let mut match_qual = [[0u8; QUALITY_TABLE_LEN]; QUALITY_TABLE_LEN];
        let mut mismatch_qual = [[0u8; QUALITY_TABLE_LEN]; QUALITY_TABLE_LEN];

        let mut q = MIN_EFFECTIVE_PHRED_INPUT;
        while usize::from(q - MIN_EFFECTIVE_PHRED_INPUT) < QUALITY_TABLE_LEN {
            error_prob[usize::from(q - MIN_EFFECTIVE_PHRED_INPUT)] = phred_to_error_prob(q);
            q = q.saturating_add(1);
        }

        let mut fwd_q = 0usize;
        while fwd_q < QUALITY_TABLE_LEN {
            let mut rev_q = 0usize;
            while rev_q < QUALITY_TABLE_LEN {
                match_qual[fwd_q][rev_q] =
                    corrected_match_quality(error_prob[fwd_q], error_prob[rev_q]);
                mismatch_qual[fwd_q][rev_q] =
                    corrected_mismatch_quality(error_prob[fwd_q], error_prob[rev_q]);
                rev_q += 1;
            }
            fwd_q += 1;
        }

        Self {
            match_qual,
            mismatch_qual,
        }
    }

    fn correct_overlap_column(
        &self,
        fwd_base: u8,
        rev_base: u8,
        fwd_qual: u8,
        rev_qual: u8,
        params: CorrectionParams,
    ) -> (u8, u8) {
        debug_assert!(is_canonical_base(fwd_base));
        debug_assert!(is_canonical_base(rev_base));

        let fwd_idx = qual_index(fwd_qual);
        let rev_idx = qual_index(rev_qual);

        if fwd_base == rev_base {
            (
                fwd_base,
                self.match_qual[fwd_idx][rev_idx].min(params.max_output_qual),
            )
        } else if fwd_qual >= rev_qual {
            (
                fwd_base,
                self.mismatch_qual[fwd_idx.max(rev_idx)][fwd_idx.min(rev_idx)]
                    .min(params.max_output_qual),
            )
        } else {
            (
                rev_base,
                self.mismatch_qual[rev_idx.max(fwd_idx)][rev_idx.min(fwd_idx)]
                    .min(params.max_output_qual),
            )
        }
    }
}

#[inline]
fn qual_index(qual: u8) -> usize {
    usize::from(
        qual.clamp(MIN_EFFECTIVE_PHRED_INPUT, MAX_EFFECTIVE_PHRED_INPUT)
            - MIN_EFFECTIVE_PHRED_INPUT,
    )
}

#[inline]
fn phred_to_error_prob(phred: u8) -> f64 {
    10_f64.powf(-f64::from(phred) / 10.0)
}

#[inline]
fn posterior_to_quality(posterior: f64) -> u8 {
    let score = (posterior.log10() * -10.0).floor();

    if !score.is_finite() || score <= 0.0 {
        return 0;
    }

    let mut quality = 0;
    while quality < MAX_CORRECTED_PHRED_OUTPUT && f64::from(quality + 1) <= score {
        quality += 1;
    }
    quality
}

#[inline]
fn corrected_match_quality(fwd_error: f64, rev_error: f64) -> u8 {
    posterior_to_quality(mismatch_error_probability(fwd_error, rev_error))
}

#[inline]
fn corrected_mismatch_quality(fwd_error: f64, rev_error: f64) -> u8 {
    posterior_to_quality(match_error_probability(fwd_error, rev_error))
}

#[inline]
fn is_canonical_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use crate::{
        merge::MergedConsensus,
        overlap::{OrientedPairSlices, OverlapBounds},
        prelude::utils::decode_fastq_quality_scores,
        test_fixtures::TupleRecord,
    };
    use proptest::prelude::*;

    fn merged_fixture<'a>(
        id: &'a str,
        seq: &[u8],
        qual: &[u8],
        fwd_source_seq: &'a [u8],
        fwd_source_qual: &[u8],
        rev_source_seq: &[u8],
        rev_source_qual: &[u8],
    ) -> (MergedConsensus, PairOverlap<'a>) {
        let left_overhang_len = seq.len().saturating_sub(fwd_source_seq.len());

        let consensus = MergedConsensus::try_new(
            id.to_string(),
            seq.to_vec(),
            decode_fastq_quality_scores(qual).into_vec(),
            left_overhang_len,
        )
        .expect("merged correction fixture should have consistent consensus lengths");

        let slices = OrientedPairSlices {
            id,
            fwd_seq: fwd_source_seq,
            fwd_qual: decode_fastq_quality_scores(fwd_source_qual),
            rev_seq_rc: rev_source_seq.to_vec().into_boxed_slice(),
            rev_qual_rev: decode_fastq_quality_scores(rev_source_qual),
        };
        let overlap = PairOverlap::from_oriented_slices(
            slices,
            OverlapBounds::new(fwd_source_seq.len(), 0, 0),
        )
        .expect("merged correction fixture should have consistent overlap slices");

        (consensus, overlap)
    }

    #[test]
    fn test_compute_corrected_score_prefers_higher_quality_on_mismatch() {
        let fwd_base = b'A';
        let rev_base = b'C';
        let fwd_qual = 35_u8;
        let rev_qual = 20_u8;

        let (base, qual) = CORRECTION_TABLES.correct_overlap_column(
            fwd_base,
            rev_base,
            fwd_qual,
            rev_qual,
            CorrectionParams::default(),
        );

        assert_eq!(base, b'A');
        assert!(qual <= 40);
    }

    #[test]
    fn test_compute_corrected_score_returns_input_base_on_match() {
        let fwd_base = b'G';
        let rev_base = b'G';
        let fwd_qual = 30_u8;
        let rev_qual = 30_u8;

        let (base, qual) = CORRECTION_TABLES.correct_overlap_column(
            fwd_base,
            rev_base,
            fwd_qual,
            rev_qual,
            CorrectionParams::default(),
        );

        assert_eq!(base, b'G');
        assert!(qual <= 40);
    }

    #[test]
    fn test_table_driven_correction_matches_scalar_oracle() {
        for fwd_qual in 0u8..=40u8 {
            for rev_qual in 0u8..=40u8 {
                let (match_base, match_qual) = CORRECTION_TABLES.correct_overlap_column(
                    b'A',
                    b'A',
                    fwd_qual,
                    rev_qual,
                    CorrectionParams::default(),
                );
                assert_eq!(match_base, b'A');
                assert_eq!(
                    match_qual,
                    scalar_corrected_quality(true, fwd_qual, rev_qual)
                );

                let (mismatch_base, mismatch_qual) = CORRECTION_TABLES.correct_overlap_column(
                    b'A',
                    b'C',
                    fwd_qual,
                    rev_qual,
                    CorrectionParams::default(),
                );
                let expected_base = if fwd_qual >= rev_qual { b'A' } else { b'C' };
                assert_eq!(mismatch_base, expected_base);
                assert_eq!(
                    mismatch_qual,
                    scalar_corrected_quality(false, fwd_qual.max(rev_qual), fwd_qual.min(rev_qual))
                );
            }
        }
    }

    #[test]
    fn test_correct_preserves_id_and_sequence() {
        let (consensus, overlap) = merged_fixture(
            "read1", b"ACGT", b"IIII", b"ACGT", b"IIII", b"ACGT", b"IIII",
        );

        let corrected = OverlapCorrector::new(CorrectionParams::default())
            .correct_merged_consensus(consensus, &overlap)
            .expect("correction should succeed for a fully consistent synthetic merged read");
        assert_eq!(corrected.id(), "read1");
        assert_eq!(corrected.sequence_bytes(), b"ACGT");
        assert_eq!(
            corrected.sequence_bytes().len(),
            corrected.quality_score_bytes().len()
        );
    }

    #[test]
    fn test_corrected_merged_into_record_roundtrip() {
        let (consensus, overlap) = merged_fixture(
            "read-merged",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
        );

        let record: TupleRecord = OverlapCorrector::new(CorrectionParams::default())
            .correct_merged_consensus(consensus, &overlap)
            .expect("correction should succeed before converting to a record")
            .into_record()
            .expect("corrected merged read should convert into a tuple record");
        assert_eq!(record.id(), "read-merged");
        assert_eq!(record.seq(), "ACGT");
        assert_eq!(record.qual(), "IIII");
    }

    #[test]
    fn test_corrected_pair_into_records_roundtrip() {
        let corrected = CorrectedOrientedPair {
            id: "read-pair".to_string(),
            fwd_seq: b"AAAA".to_vec(),
            fwd_quality_scores: vec![40; 4],
            rev_seq_rc: b"AAAA".to_vec(),
            rev_quality_scores_rc: vec![41; 4],
            overlap_bounds: OverlapBounds::new(4, 0, 0),
        };

        let (left, right) = OwnedReadPair::try_from(corrected)
            .expect("corrected pair should convert into owned reads")
            .into_reads();
        let left = TupleRecord::try_from(left)
            .expect("left corrected read should convert into tuple record");
        let right = TupleRecord::try_from(right)
            .expect("right corrected read should convert into tuple record");
        assert_eq!(left.id(), "read-pair");
        assert_eq!(left.seq(), "AAAA");
        assert_eq!(left.qual(), "IIII");
        assert_eq!(right.id(), "read-pair");
        assert_eq!(right.seq(), "TTTT");
        assert_eq!(right.qual(), "JJJJ");
    }

    #[test]
    fn test_corrected_qualities_match_consensus_len_with_overhangs() {
        let (consensus, overlap) = merged_fixture(
            "read1",
            b"TTTTACGT",
            b"IIIIIIII",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
        );

        let corrected = OverlapCorrector::new(CorrectionParams::default())
            .correct_merged_consensus(consensus, &overlap)
            .expect("correction should not error for overhang-quality regression fixture");
        assert_eq!(
            corrected.sequence_bytes().len(),
            corrected.quality_score_bytes().len()
        );
    }

    #[test]
    fn test_correct_preserves_non_overlap_qualities_in_merged_output() {
        let consensus = MergedConsensus::try_new(
            "read-overhangs".to_string(),
            b"TTTTACGTGG".to_vec(),
            decode_fastq_quality_scores(b"JKLMIIIIWX").into_vec(),
            4,
        )
        .expect("merged overhang-preservation fixture should have consistent layout");
        let overlap = PairOverlap::from_oriented_slices(
            OrientedPairSlices {
                id: "read-overhangs",
                fwd_seq: b"ACGT",
                fwd_qual: decode_fastq_quality_scores(b"IIII"),
                rev_seq_rc: b"ACGT".to_vec().into_boxed_slice(),
                rev_qual_rev: decode_fastq_quality_scores(b"IIII"),
            },
            OverlapBounds::new(4, 0, 0),
        )
        .expect("merged overhang-preservation fixture should have consistent overlap");

        let corrected = OverlapCorrector::new(CorrectionParams::default())
            .correct_merged_consensus(consensus, &overlap)
            .expect("correction should succeed for merged overhang-preservation fixture");

        assert_eq!(&corrected.quality_score_bytes()[..4], [41, 42, 43, 44]);
        assert_eq!(&corrected.quality_score_bytes()[8..], [54, 55]);
    }

    #[test]
    fn test_max_output_qual_caps_correction_scores() {
        let (_, uncapped) = CORRECTION_TABLES.correct_overlap_column(
            b'A',
            b'A',
            40,
            40,
            CorrectionParams::default(),
        );
        assert!(uncapped > 10);

        let window = CorrectionWindow::new(b"A", &[40], b"A", &[40]);
        let mut seq = [b'A'];
        let mut qual = [0u8];

        OverlapCorrector::new(CorrectionParams::default().with_max_output_qual(10))
            .correct_merged_overlap(&window, &mut seq, &mut qual);

        assert_eq!(qual[0], 10);
    }

    #[test]
    fn test_quality_only_preserves_forward_base_choice_on_mismatch() {
        let overlap = PairOverlap::from_oriented_slices(
            OrientedPairSlices {
                id: "read1",
                fwd_seq: b"A",
                fwd_qual: [0].into(),
                rev_seq_rc: b"G".as_slice().into(),
                rev_qual_rev: [40].into(),
            },
            OverlapBounds::new(1, 0, 0),
        )
        .expect("single-base overlap fixture should be valid");

        let corrected = OverlapCorrector::new(CorrectionParams::default().quality_only())
            .correct_pair_overlap(&overlap)
            .expect("correcting from pair-overlap slices should succeed");
        let (left, right) = OwnedReadPair::try_from(corrected)
            .expect("corrected pair should convert to owned reads")
            .into_reads();

        assert_eq!(left.sequence_bytes(), b"A");
        assert_eq!(right.sequence_bytes(), b"T");
    }

    proptest! {
        #[test]
        fn proptest_compute_corrected_score_respects_basic_kernel_invariants(
            fwd_base in prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            rev_base in prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            fwd_qual in 0u8..=40u8,
            rev_qual in 0u8..=40u8,
        ) {
            let (chosen_base, corrected_qual) = CORRECTION_TABLES.correct_overlap_column(
                fwd_base,
                rev_base,
                fwd_qual,
                rev_qual,
                CorrectionParams::default(),
            );

            prop_assert!(chosen_base == fwd_base || chosen_base == rev_base);
            prop_assert!(corrected_qual <= 40);

            if fwd_base == rev_base {
                prop_assert_eq!(chosen_base, fwd_base);
            } else if fwd_qual >= rev_qual {
                prop_assert_eq!(chosen_base, fwd_base);
            } else {
                prop_assert_eq!(chosen_base, rev_base);
            }
        }
    }

    fn scalar_corrected_quality(is_match: bool, chosen_qual: u8, other_qual: u8) -> u8 {
        let chosen_error = phred_to_error_prob(
            chosen_qual.clamp(MIN_EFFECTIVE_PHRED_INPUT, MAX_EFFECTIVE_PHRED_INPUT),
        );
        let other_error = phred_to_error_prob(
            other_qual.clamp(MIN_EFFECTIVE_PHRED_INPUT, MAX_EFFECTIVE_PHRED_INPUT),
        );
        let posterior = if is_match {
            mismatch_error_probability(chosen_error, other_error)
        } else {
            match_error_probability(chosen_error, other_error)
        };

        posterior_to_quality(posterior)
    }
}
