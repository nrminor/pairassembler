use std::{borrow::Cow, fmt::Display};

use crate::{
    PairOverlap, ReadPair, Result,
    assembler::{
        FromRecordParts, IntoOwnedPairRecordParts, IntoOwnedRecordParts, IntoRecordConversion,
        IntoRecordsConversion,
    },
    errors::CorrectionError::ConsensusLengthMismatch,
    merge::MergedRead,
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Configuration for correction behavior.
///
/// This remains a placeholder until the correction module is redesigned; it exists so the
/// assembler API can reserve a stable place for future correction controls.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CorrectionParams {}

/// Aligned overlap-local input window shared by correction kernels.
#[derive(Debug, Clone)]
pub struct CorrectionWindow<'a> {
    fwd_seq: Cow<'a, [u8]>,
    fwd_qual: Cow<'a, [u8]>,
    rev_seq: Cow<'a, [u8]>,
    rev_qual: Cow<'a, [u8]>,
}

/// Corrected consensus record emitted after applying overlap-based quality correction.
#[derive(Debug)]
pub struct CorrectedMergedRead {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Corrected paired reads emitted by pair-preserving correction flows.
#[derive(Debug)]
pub struct CorrectedReadPair {
    id: String,
    fwd_seq: Vec<u8>,
    fwd_qual: Vec<u8>,
    rev_seq: Vec<u8>,
    rev_qual: Vec<u8>,
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
    pub fn fwd_quality_bytes(&self) -> &[u8] {
        self.fwd_qual.as_slice()
    }

    #[must_use]
    pub fn rev_sequence_bytes(&self) -> &[u8] {
        self.rev_seq.as_slice()
    }

    #[must_use]
    pub fn rev_quality_bytes(&self) -> &[u8] {
        self.rev_qual.as_slice()
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
        fwd_seq: Cow<'a, [u8]>,
        fwd_qual: Cow<'a, [u8]>,
        rev_seq: Cow<'a, [u8]>,
        rev_qual: Cow<'a, [u8]>,
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
    pub fn len(&self) -> usize {
        self.fwd_seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn forward_sequence(&self) -> &[u8] {
        &self.fwd_seq
    }

    #[must_use]
    pub fn forward_qualities(&self) -> &[u8] {
        &self.fwd_qual
    }

    #[must_use]
    pub fn reverse_sequence(&self) -> &[u8] {
        &self.rev_seq
    }

    #[must_use]
    pub fn reverse_qualities(&self) -> &[u8] {
        &self.rev_qual
    }
}

impl CorrectedMergedRead {
    /// Construct a corrected merged read from checked record parts.
    ///
    /// # Errors
    ///
    /// Returns an error if sequence and quality lengths differ.
    pub(crate) fn try_new(id: String, seq: Vec<u8>, qual: Vec<u8>) -> Result<Self> {
        if seq.len() != qual.len() {
            return Err(ConsensusLengthMismatch {
                seq_len: seq.len(),
                qual_len: qual.len(),
            }
            .into());
        }

        Ok(Self { id, seq, qual })
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
    pub fn quality_bytes(&self) -> &[u8] {
        self.qual.as_slice()
    }

    #[must_use]
    pub fn qualities_owned(self) -> Vec<u8> {
        self.qual
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
        (self.id, self.seq, self.qual)
    }
}

impl MergedRead {
    /// Apply quality score correction across the merged overlap.
    ///
    /// # Errors
    ///
    /// Returns an error if the corrected quality vector cannot be reconciled with the existing
    /// merged consensus layout.
    pub fn correct(self) -> Result<CorrectedMergedRead> {
        // Pull out the ID and the sequence from prior to correction, as we'll be recycling these.
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

        // Run correction on the quality scores, for which we'll use a handy parallel iterator from rayon
        let window = CorrectionWindow::new(
            Cow::Owned(fwd_source_seq),
            Cow::Owned(fwd_source_qual),
            Cow::Owned(rev_source_seq),
            Cow::Owned(rev_source_qual),
        );
        let overlap_corrected_quals = correct_window(&window)
            .into_iter()
            .map(|(_, qual)| qual)
            .collect::<Vec<_>>();

        // Preserve non-overlap qualities and overwrite only the merged overlap window.
        let overlap_end = left_overhang_len + overlap_corrected_quals.len();
        let mut corrected_quals = consensus_qual;
        corrected_quals[left_overhang_len..overlap_end].copy_from_slice(&overlap_corrected_quals);

        CorrectedMergedRead::try_new(id, seq, corrected_quals)
    }
}

impl ReadPair<'_> {
    pub(crate) fn correct_from_overlap(&self, overlap: &PairOverlap) -> CorrectedReadPair {
        let mut fwd_seq = self.fwd_sequence_bytes().to_vec();
        let mut rev_seq = self.rev_sequence_bytes().to_vec();
        let mut fwd_qual = self
            .fwd_quality_bytes()
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect::<Vec<_>>();
        let mut rev_qual = self
            .rev_quality_bytes()
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect::<Vec<_>>();

        let fwd_qual_ascii = self.fwd_quality_bytes();

        let window = CorrectionWindow::new(
            Cow::Borrowed(overlap.forward_sequence()),
            Cow::Borrowed(overlap.forward_qualities()),
            Cow::Borrowed(overlap.reverse_sequence()),
            Cow::Borrowed(overlap.reverse_qualities()),
        );
        let corrected_overlap = correct_window(&window);

        for (i, (chosen_base_rc, corrected_q)) in corrected_overlap.into_iter().enumerate() {
            let fwd_idx = overlap.forward_start_offset() + i;
            let rev_rc_idx = overlap.reverse_start_offset() + i;
            let rev_idx = self.rev_mate.len() - 1 - rev_rc_idx;

            fwd_seq[fwd_idx] = chosen_base_rc;
            fwd_qual[fwd_idx] = corrected_q;
            rev_seq[rev_idx] = complement_base(chosen_base_rc);
            rev_qual[rev_idx] = corrected_q;
        }

        CorrectedReadPair {
            id: self.fwd_id().to_string(),
            fwd_seq,
            fwd_qual,
            rev_seq,
            rev_qual,
        }
    }
}

fn correct_window(window: &CorrectionWindow<'_>) -> Vec<(u8, u8)> {
    (0..window.len())
        .into_par_iter()
        .map(|idx| {
            let fwd_base = window.forward_sequence()[idx];
            let rev_base = window.reverse_sequence()[idx];
            let fwd_qual = window.forward_qualities()[idx];
            let rev_qual = window.reverse_qualities()[idx];

            let base_overlap = BaseOverlap {
                fwd_base: &fwd_base,
                rev_base: &rev_base,
                fwd_qual: &fwd_qual,
                rev_qual: &rev_qual,
            };

            let (base, qual) = base_overlap.compute_corrected_score();
            (*base, qual)
        })
        .collect::<Vec<_>>()
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

#[derive(Debug, Clone)]
pub struct BaseOverlap<'overlap> {
    fwd_base: &'overlap u8,
    rev_base: &'overlap u8,
    fwd_qual: &'overlap u8,
    rev_qual: &'overlap u8,
}

impl<'overlap> BaseOverlap<'overlap> {
    pub fn new(
        fwd_base: &'overlap u8,
        rev_base: &'overlap u8,
        fwd_qual: &'overlap u8,
        rev_qual: &'overlap u8,
    ) -> Self {
        Self {
            fwd_base,
            rev_base,
            fwd_qual,
            rev_qual,
        }
    }

    // TODO: This may need to be modified to support correction of unmerged reads
    pub fn compute_corrected_score(&self) -> (&'overlap u8, u8) {
        // run some checks if in debug mode before proceeding
        debug_assert!(
            self.fwd_qual.saturating_sub(33) <= 60 && self.rev_qual.saturating_sub(33) <= 60,
            "Unusually high quality scores detected"
        );
        debug_assert!(
            matches!(*self.fwd_base, b'A' | b'C' | b'G' | b'T'),
            "Unexpected base in forward read: {}",
            *self.fwd_base as char
        );
        debug_assert!(
            matches!(*self.rev_base, b'A' | b'C' | b'G' | b'T'),
            "Unexpected base in reverse read: {}",
            *self.rev_base as char
        );

        // run some casts for more precision and convert the Phred score into an error likelihood
        let fwd_qual = f64::from(self.fwd_qual.saturating_sub(33));
        let rev_qual = f64::from(self.rev_qual.saturating_sub(33));
        let fwd_error = 10_f64.powf(-fwd_qual / 10.0);
        let rev_error = 10_f64.powf(-rev_qual / 10.0);

        if self.fwd_base == self.rev_base {
            let status = Match {
                fwd_error: &fwd_error,
                rev_error: &rev_error,
            };
            let score = status.compute_score();
            (self.fwd_base, score)
        } else {
            let status = Mismatch {
                fwd_error: &fwd_error,
                rev_error: &rev_error,
            };
            let score = status.compute_score();
            if self.fwd_qual >= self.rev_qual {
                (self.fwd_base, score)
            } else {
                (self.rev_base, score)
            }
        }
    }
}

pub enum MatchStatus<'err_prob> {
    Match {
        fwd_error: &'err_prob f64,
        rev_error: &'err_prob f64,
    },
    Mismatch {
        fwd_error: &'err_prob f64,
        rev_error: &'err_prob f64,
    },
}
pub use MatchStatus::*;

impl MatchStatus<'_> {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn compute_score(self) -> u8 {
        let posterior = match self {
            Match {
                fwd_error,
                rev_error,
            } => mismatch_error_probability(*fwd_error, *rev_error),
            Mismatch {
                fwd_error,
                rev_error,
            } => match_error_probability(*fwd_error, *rev_error),
        };

        // compute the integer quality score
        let score = (posterior.log10() * -10.0).floor();

        if score > 40.0 { 40_u8 } else { score as u8 }
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

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use crate::merge::MergeProvenance;
    use crate::test_fixtures::TupleRecord;

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
            fwd_source_qual.to_vec(),
            rev_source_seq.to_vec(),
            rev_source_qual.to_vec(),
        )
        .expect("merged correction fixture should have consistent provenance lengths");

        let left_overhang_len = seq.len().saturating_sub(fwd_source_seq.len());

        MergedRead::try_new(
            id.to_string(),
            seq.to_vec(),
            qual.to_vec(),
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

        let overlap = BaseOverlap::new(&fwd_base, &rev_base, &fwd_qual, &rev_qual);
        let (base, qual) = overlap.compute_corrected_score();

        assert_eq!(*base, b'A');
        assert!(qual <= 40);
    }

    #[test]
    fn test_compute_corrected_score_returns_input_base_on_match() {
        let fwd_base = b'G';
        let rev_base = b'G';
        let fwd_qual = 30_u8;
        let rev_qual = 30_u8;

        let overlap = BaseOverlap::new(&fwd_base, &rev_base, &fwd_qual, &rev_qual);
        let (base, qual) = overlap.compute_corrected_score();

        assert_eq!(*base, b'G');
        assert!(qual <= 40);
    }

    #[test]
    fn test_correct_preserves_id_and_sequence() {
        let uncorrected = merged_fixture(
            "read1", b"ACGT", b"IIII", b"ACGT", b"IIII", b"ACGT", b"IIII",
        );

        let corrected = uncorrected
            .correct()
            .expect("correction should succeed for a fully consistent synthetic merged read");
        assert_eq!(corrected.id(), "read1");
        assert_eq!(corrected.sequence_bytes(), b"ACGT");
        assert_eq!(
            corrected.sequence_bytes().len(),
            corrected.quality_bytes().len()
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
            .correct()
            .expect("correction should succeed before converting to a record")
            .into_record()
            .expect("corrected merged read should convert into a tuple record");
        assert_eq!(record.id(), "read-merged");
        assert_eq!(record.seq(), "ACGT");
        assert_eq!(record.qual(), "((((");
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
            .correct()
            .expect("correction should not error for overhang-quality regression fixture");
        assert_eq!(
            corrected.sequence_bytes().len(),
            corrected.quality_bytes().len()
        );
    }
}
