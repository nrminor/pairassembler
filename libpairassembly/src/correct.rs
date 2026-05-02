use std::{fmt::Display, sync::LazyLock};

use crate::{
    PairOverlap, ReadPair, Result,
    assembler::{
        FromRecordParts, IntoOwnedPairRecordParts, IntoOwnedRecordParts, IntoRecordConversion,
        IntoRecordsConversion,
    },
    errors::CorrectionError::ConsensusLengthMismatch,
    merge::{MergedConsensus, MergedRead},
    prelude::utils::{encode_fastq_quality_scores_in_place, fastq_ascii_to_phred},
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

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
/// This surface is intentionally small for now and marked non-exhaustive because the correction
/// policy is still evolving.
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
pub struct CorrectionWindow<'a> {
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

/// Corrected paired reads emitted by pair-preserving correction flows.
#[derive(Debug, Clone)]
pub struct CorrectedReadPair {
    id: String,
    fwd_seq: Vec<u8>,
    fwd_qual: Vec<u8>,
    rev_seq: Vec<u8>,
    rev_qual: Vec<u8>,
}

struct MergedCorrectionInput<'a> {
    id: String,
    seq: Vec<u8>,
    consensus_qual: Vec<u8>,
    left_overhang_len: usize,
    fwd_source_seq: Vec<u8>,
    fwd_source_qual: Vec<u8>,
    rev_source_seq: Vec<u8>,
    rev_source_qual: Vec<u8>,
    window: CorrectionWindow<'a>,
}

struct PairCorrectionInput<'a> {
    id: String,
    fwd_seq: Vec<u8>,
    fwd_qual: Vec<u8>,
    rev_seq: Vec<u8>,
    rev_qual: Vec<u8>,
    overlap: &'a PairOverlap<'a>,
    window: CorrectionWindow<'a>,
}

impl CorrectedReadPair {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn fwd_sequence_bytes(&self) -> &[u8] {
        self.fwd_seq.as_slice()
    }

    #[must_use]
    pub fn fwd_sequence(&self) -> &str {
        std::str::from_utf8(self.fwd_sequence_bytes())
            .expect("corrected forward sequence should remain valid ASCII DNA")
    }

    #[must_use]
    pub fn fwd_quality_bytes(&self) -> &[u8] {
        self.fwd_qual.as_slice()
    }

    #[must_use]
    pub fn fwd_quality_scores(&self) -> &str {
        std::str::from_utf8(self.fwd_quality_bytes())
            .expect("corrected forward qualities should remain valid ASCII")
    }

    #[must_use]
    pub fn rev_sequence_bytes(&self) -> &[u8] {
        self.rev_seq.as_slice()
    }

    #[must_use]
    pub fn rev_sequence(&self) -> &str {
        std::str::from_utf8(self.rev_sequence_bytes())
            .expect("corrected reverse sequence should remain valid ASCII DNA")
    }

    #[must_use]
    pub fn rev_quality_bytes(&self) -> &[u8] {
        self.rev_qual.as_slice()
    }

    #[must_use]
    pub fn rev_quality_scores(&self) -> &str {
        std::str::from_utf8(self.rev_quality_bytes())
            .expect("corrected reverse qualities should remain valid ASCII")
    }

    /// Convert corrected paired output into two user record values.
    ///
    /// This is a boundary conversion API. Identity-shaped `into_*` methods are
    /// intentionally omitted to keep `into_*` naming reserved for meaningful
    /// representation changes.
    ///
    /// # Errors
    ///
    /// Returns an error if either target record type cannot be constructed
    /// from the corrected parts.
    pub fn into_records<T>(self) -> Result<(T, T)>
    where
        T: FromRecordParts,
        T::Error: Display,
    {
        IntoRecordsConversion::into_records(self)
    }
}

impl<'a> CorrectionWindow<'a> {
    #[must_use]
    pub(crate) fn new(
        fwd_seq: &'a [u8],
        fwd_qual: &'a [u8],
        rev_seq: &'a [u8],
        rev_qual: &'a [u8],
    ) -> Self {
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
    pub(crate) fn from_overlap<'pair>(overlap: &'a PairOverlap<'pair>) -> Self {
        Self::new(
            overlap.forward_sequence(),
            overlap.forward_qualities(),
            overlap.reverse_sequence(),
            overlap.reverse_qualities(),
        )
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.fwd_seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn forward_sequence(&self) -> &[u8] {
        self.fwd_seq
    }

    #[must_use]
    pub fn forward_qualities(&self) -> &[u8] {
        self.fwd_qual
    }

    #[must_use]
    pub fn reverse_sequence(&self) -> &[u8] {
        self.rev_seq
    }

    #[must_use]
    pub fn reverse_qualities(&self) -> &[u8] {
        self.rev_qual
    }

    fn correct_merged_in_place(
        &self,
        tables: &CorrectionTables,
        params: CorrectionParams,
        seq_overlap: &mut [u8],
        qual_overlap: &mut [u8],
    ) {
        debug_assert_eq!(seq_overlap.len(), self.len());
        debug_assert_eq!(qual_overlap.len(), self.len());

        for idx in 0..self.len() {
            let (corrected_base, corrected_qual) = tables.correct_overlap_column(
                self.forward_sequence()[idx],
                self.reverse_sequence()[idx],
                self.forward_qualities()[idx],
                self.reverse_qualities()[idx],
                params,
            );

            if !params.quality_only {
                seq_overlap[idx] = corrected_base;
            }
            qual_overlap[idx] = corrected_qual;
        }
    }

    fn correct_pair_in_place(
        &self,
        tables: &CorrectionTables,
        params: CorrectionParams,
        fwd_seq_overlap: &mut [u8],
        fwd_qual_overlap: &mut [u8],
        rev_seq_overlap_raw: &mut [u8],
        rev_qual_overlap_raw: &mut [u8],
    ) {
        debug_assert_eq!(fwd_seq_overlap.len(), self.len());
        debug_assert_eq!(fwd_qual_overlap.len(), self.len());
        debug_assert_eq!(rev_seq_overlap_raw.len(), self.len());
        debug_assert_eq!(rev_qual_overlap_raw.len(), self.len());

        for idx in 0..self.len() {
            let (corrected_base, corrected_qual) = tables.correct_overlap_column(
                self.forward_sequence()[idx],
                self.reverse_sequence()[idx],
                self.forward_qualities()[idx],
                self.reverse_qualities()[idx],
                params,
            );

            let write_base = if params.quality_only {
                self.forward_sequence()[idx]
            } else {
                corrected_base
            };

            fwd_seq_overlap[idx] = write_base;
            fwd_qual_overlap[idx] = corrected_qual;

            let rev_raw_idx = self.len() - 1 - idx;
            rev_seq_overlap_raw[rev_raw_idx] = complement_base(write_base);
            rev_qual_overlap_raw[rev_raw_idx] = corrected_qual;
        }
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
        T: FromRecordParts,
        T::Error: Display,
    {
        IntoRecordConversion::into_record(self)
    }

    /// Correct a score-space merged consensus using retained overlap evidence.
    ///
    /// # Errors
    ///
    /// Returns an error if the corrected quality vector cannot be reconciled with the existing
    /// merged consensus layout.
    pub(crate) fn correct_consensus_with(
        consensus: MergedConsensus,
        overlap: &PairOverlap<'_>,
        params: CorrectionParams,
    ) -> Result<Self> {
        let window = CorrectionWindow::from_overlap(overlap);
        let overlap_start = consensus.left_overhang_len();
        let overlap_end = overlap_start + window.len();
        let (id, mut corrected_seq, mut corrected_quals) = consensus.into_score_parts();

        window.correct_merged_in_place(
            &CORRECTION_TABLES,
            params,
            &mut corrected_seq[overlap_start..overlap_end],
            &mut corrected_quals[overlap_start..overlap_end],
        );

        Self::try_new(id, corrected_seq, corrected_quals)
    }
}

impl IntoOwnedPairRecordParts for CorrectedReadPair {
    fn into_owned_pair_record_parts(self) -> (String, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
        (
            self.id,
            self.fwd_seq,
            self.fwd_qual,
            self.rev_seq,
            self.rev_qual,
        )
    }
}

impl IntoOwnedRecordParts for CorrectedMergedRead {
    fn into_owned_record_parts(self) -> (String, Vec<u8>, Vec<u8>) {
        let mut quality_ascii = self.quality_scores;
        encode_fastq_quality_scores_in_place(&mut quality_ascii);

        (self.id, self.seq, quality_ascii)
    }
}

impl MergedRead {
    /// Apply correction across the merged overlap with explicit correction params.
    ///
    /// # Errors
    ///
    /// Returns an error if the corrected quality vector cannot be reconciled with the existing
    /// merged consensus layout.
    pub fn correct_with(self, params: CorrectionParams) -> Result<CorrectedMergedRead> {
        let (
            id,
            seq,
            consensus_qual,
            left_overhang_len,
            fwd_source_seq,
            fwd_source_qual,
            rev_source_seq,
            rev_source_qual,
        ) = self.into_correction_parts();

        let window = CorrectionWindow::new(
            &fwd_source_seq,
            &fwd_source_qual,
            &rev_source_seq,
            &rev_source_qual,
        );
        let overlap_end = left_overhang_len + window.len();
        let mut corrected_quals = consensus_qual;
        let mut corrected_seq = seq;

        window.correct_merged_in_place(
            &CORRECTION_TABLES,
            params,
            &mut corrected_seq[left_overhang_len..overlap_end],
            &mut corrected_quals[left_overhang_len..overlap_end],
        );
        CorrectedMergedRead::try_new(id, corrected_seq, corrected_quals)
    }
}

impl ReadPair<'_> {
    pub(crate) fn correct_from_overlap_with(
        &self,
        overlap: &PairOverlap,
        params: CorrectionParams,
    ) -> CorrectedReadPair {
        let mut fwd_seq = self.fwd_sequence_bytes().to_vec();
        let mut rev_seq = self.rev_sequence_bytes().to_vec();
        let mut fwd_qual = self
            .fwd_quality_bytes()
            .iter()
            .map(|q| fastq_ascii_to_phred(*q))
            .collect::<Vec<_>>();
        let mut rev_qual = self
            .rev_quality_bytes()
            .iter()
            .map(|q| fastq_ascii_to_phred(*q))
            .collect::<Vec<_>>();

        let window = CorrectionWindow::from_overlap(overlap);
        let fwd_start = overlap.forward_start_offset();
        let fwd_end = fwd_start + window.len();
        let rev_raw_start = rev_seq.len() - (overlap.reverse_start_offset() + window.len());
        let rev_raw_end = rev_seq.len() - overlap.reverse_start_offset();

        window.correct_pair_in_place(
            &CORRECTION_TABLES,
            params,
            &mut fwd_seq[fwd_start..fwd_end],
            &mut fwd_qual[fwd_start..fwd_end],
            &mut rev_seq[rev_raw_start..rev_raw_end],
            &mut rev_qual[rev_raw_start..rev_raw_end],
        );
        encode_fastq_quality_scores_in_place(&mut fwd_qual);
        encode_fastq_quality_scores_in_place(&mut rev_qual);

        CorrectedReadPair {
            id: self.fwd_id().to_string(),
            fwd_seq,
            fwd_qual,
            rev_seq,
            rev_qual,
        }
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
    error_prob: [f64; QUALITY_TABLE_LEN],
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
            error_prob,
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

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
#[inline]
fn posterior_to_quality(posterior: f64) -> u8 {
    let score = (posterior.log10() * -10.0).floor();
    if score > f64::from(MAX_CORRECTED_PHRED_OUTPUT) {
        MAX_CORRECTED_PHRED_OUTPUT
    } else {
        score as u8
    }
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
    use crate::SequenceRead;
    use crate::merge::MergeProvenance;
    use crate::overlap::{OverlapBounds, PreparedPair};
    use crate::prelude::utils::decode_fastq_quality_scores;
    use crate::test_fixtures::TupleRecord;
    use proptest::prelude::*;

    fn merged_fixture(
        id: &str,
        seq: &[u8],
        qual: &[u8],
        fwd_source_seq: &[u8],
        fwd_source_qual: &[u8],
        rev_source_seq: &[u8],
        rev_source_qual: &[u8],
    ) -> MergedRead {
        let provenance = MergeProvenance::try_new(
            fwd_source_seq.len(),
            fwd_source_seq.to_vec(),
            decode_fastq_quality_scores(fwd_source_qual).into_vec(),
            rev_source_seq.to_vec(),
            decode_fastq_quality_scores(rev_source_qual).into_vec(),
        )
        .expect("merged correction fixture should have consistent provenance lengths");

        let left_overhang_len = seq.len().saturating_sub(fwd_source_seq.len());

        MergedRead::try_new(
            id.to_string(),
            seq.to_vec(),
            decode_fastq_quality_scores(qual).into_vec(),
            left_overhang_len,
            provenance,
        )
        .expect("merged correction fixture should have consistent consensus lengths")
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
        let uncorrected = merged_fixture(
            "read1", b"ACGT", b"IIII", b"ACGT", b"IIII", b"ACGT", b"IIII",
        );

        let corrected = uncorrected
            .correct_with(CorrectionParams::default())
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
        let uncorrected = merged_fixture(
            "read-merged",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
        );

        let record: TupleRecord = uncorrected
            .correct_with(CorrectionParams::default())
            .expect("correction should succeed before converting to a record")
            .into_record()
            .expect("corrected merged read should convert into a tuple record");
        assert_eq!(record.id(), "read-merged");
        assert_eq!(record.seq(), "ACGT");
        assert_eq!(record.qual(), "IIII");
    }

    #[test]
    fn test_corrected_pair_into_records_roundtrip() {
        let corrected = CorrectedReadPair {
            id: "read-pair".to_string(),
            fwd_seq: b"AAAA".to_vec(),
            fwd_qual: b"IIII".to_vec(),
            rev_seq: b"TTTT".to_vec(),
            rev_qual: b"JJJJ".to_vec(),
        };

        let (left, right): (TupleRecord, TupleRecord) = corrected
            .into_records()
            .expect("corrected pair should convert into two tuple records");
        assert_eq!(left.id(), "read-pair");
        assert_eq!(left.seq(), "AAAA");
        assert_eq!(left.qual(), "IIII");
        assert_eq!(right.id(), "read-pair");
        assert_eq!(right.seq(), "TTTT");
        assert_eq!(right.qual(), "JJJJ");
    }

    #[test]
    fn test_corrected_qualities_match_consensus_len_with_overhangs() {
        let uncorrected = merged_fixture(
            "read1",
            b"TTTTACGT",
            b"IIIIIIII",
            b"ACGT",
            b"IIII",
            b"ACGT",
            b"IIII",
        );

        let corrected = uncorrected
            .correct_with(CorrectionParams::default())
            .expect("correction should not error for overhang-quality regression fixture");
        assert_eq!(
            corrected.sequence_bytes().len(),
            corrected.quality_score_bytes().len()
        );
    }

    #[test]
    fn test_correct_preserves_non_overlap_qualities_in_merged_output() {
        let provenance = MergeProvenance::try_new(
            4,
            b"ACGT".to_vec(),
            decode_fastq_quality_scores(b"IIII").into_vec(),
            b"ACGT".to_vec(),
            decode_fastq_quality_scores(b"IIII").into_vec(),
        )
        .expect("merged overhang-preservation fixture should have consistent provenance");
        let uncorrected = MergedRead::try_new(
            "read-overhangs".to_string(),
            b"TTTTACGTGG".to_vec(),
            decode_fastq_quality_scores(b"JKLMIIIIWX").into_vec(),
            4,
            provenance,
        )
        .expect("merged overhang-preservation fixture should have consistent layout");

        let corrected = uncorrected
            .correct_with(CorrectionParams::default())
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

        window.correct_merged_in_place(
            &CORRECTION_TABLES,
            CorrectionParams::default().with_max_output_qual(10),
            &mut seq,
            &mut qual,
        );

        assert_eq!(qual[0], 10);
    }

    #[test]
    fn test_quality_only_preserves_forward_base_choice_on_mismatch() {
        let pair = ReadPair::from(
            SequenceRead::new("read1", "A", "!"),
            SequenceRead::new("read1", "G", "I"),
        )
        .expect("single-base mismatch fixture should pair cleanly");
        let overlap = PairOverlap::from_prepared(
            PreparedPair {
                id: "read1",
                fwd_seq: b"A",
                fwd_qual: [0].into(),
                rev_seq_rc: b"G".as_slice().into(),
                rev_qual_rev: [40].into(),
            },
            OverlapBounds::new(1, 0, 0),
        )
        .expect("single-base overlap fixture should be valid");

        let corrected =
            pair.correct_from_overlap_with(&overlap, CorrectionParams::default().quality_only());

        assert_eq!(corrected.fwd_sequence_bytes(), b"A");
        assert_eq!(corrected.rev_sequence_bytes(), b"T");
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
