use std::fmt::Display;

use crate::{
    PairOverlap, ReadPair, Result,
    assembler::{
        FromRecordParts, IntoOwnedPairRecordParts, IntoOwnedRecordParts, IntoRecordConversion,
        IntoRecordsConversion,
    },
    merge::UncorrectedMergedRead,
};
use itertools::izip;
use rayon::iter::{ParallelBridge, ParallelIterator};

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CorrectionParams {}

// TODO: this should just implement SeqRecord, no?
#[derive(Debug)]
pub struct CorrectedMergedRead {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

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

impl CorrectedMergedRead {
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

impl UncorrectedMergedRead {
    /// Apply quality score correction across the merged overlap.
    ///
    /// # Errors
    ///
    /// This currently returns `Result` for API consistency with the broader pipeline.
    /// No runtime error paths are currently produced in this implementation.
    pub fn correct(self) -> Result<CorrectedMergedRead> {
        // Pull out the ID and the sequence from prior to correction, as we'll be recycling these.
        let (id, seq, fwd_source_seq, fwd_source_qual, rev_source_seq, rev_source_qual) =
            self.into_correction_parts();

        // Run correction on the quality scores, for which we'll use a handy parallel iterator from rayon
        let corrected_quals = izip!(
                fwd_source_seq,
                rev_source_seq,
                fwd_source_qual,
                rev_source_qual,
            )
            .par_bridge() // rust is seriously magic sometimes
            .map(|(fwd_base, rev_base, fwd_qual, rev_qual)| {
                // fill the necessary information for this vertical slice of the overlap into a
                // structure for score correction
                let base_overlap = BaseOverlap {
                    fwd_base: &fwd_base,
                    rev_base: &rev_base,
                    fwd_qual: &fwd_qual,
                    rev_qual: &rev_qual,
                };

                // run the score correction and return the corrected quality score
                let (_, qual) = base_overlap.compute_corrected_score();
                qual
            })
            .collect::<Vec<_>>();

        let new_read = CorrectedMergedRead {
            id,
            seq,
            qual: corrected_quals,
        };

        Ok(new_read)
    }
}

impl ReadPair<'_> {
    pub(crate) fn correct_from_overlap(&self, overlap: &PairOverlap) -> Result<CorrectedReadPair> {
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

        for i in 0..overlap.len() {
            let fwd_idx = overlap.forward_start_offset() + i;
            let rev_rc_idx = overlap.reverse_start_offset() + i;
            let rev_idx = self.rev_mate.len() - 1 - rev_rc_idx;

            let fwd_base = fwd_seq[fwd_idx];
            let rev_base = overlap.reverse_sequence()[i];
            let fwd_q = fwd_qual_ascii[fwd_idx];
            let rev_q = overlap.reverse_qualities()[i];
            let base_overlap = BaseOverlap::new(&fwd_base, &rev_base, &fwd_q, &rev_q);

            let (chosen_base_rc, corrected_q) = base_overlap.compute_corrected_score();
            fwd_seq[fwd_idx] = *chosen_base_rc;
            fwd_qual[fwd_idx] = corrected_q;
            rev_seq[rev_idx] = complement_base(*chosen_base_rc);
            rev_qual[rev_idx] = corrected_q;
        }

        Ok(CorrectedReadPair {
            id: self.fwd_id().to_string(),
            fwd_seq,
            fwd_qual,
            rev_seq,
            rev_qual,
        })
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
    use crate::{merge::UncorrectedMergedRead, test_fixtures::TupleRecord};

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
        let uncorrected = UncorrectedMergedRead::from_parts(
            "read1".to_string(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
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
        let uncorrected = UncorrectedMergedRead::from_parts(
            "read-merged".to_string(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
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
    #[ignore = "Known issue: correction only emits overlap-length qualities"]
    fn test_corrected_qualities_match_consensus_len_with_overhangs() {
        let uncorrected = UncorrectedMergedRead::from_parts(
            "read1".to_string(),
            b"TTTTACGT".to_vec(),
            b"IIIIIIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
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
