use crate::{Result, merge::UncorrectedMergedRead};
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

impl CorrectedMergedRead {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn sequence(&self) -> &[u8] {
        self.seq.as_slice()
    }

    pub fn sequence_owned(self) -> Vec<u8> {
        self.seq
    }

    pub fn qualities(&self) -> &[u8] {
        self.qual.as_slice()
    }

    pub fn qualities_owned(self) -> Vec<u8> {
        self.qual
    }
}

impl UncorrectedMergedRead<'_> {
    pub fn correct(self) -> Result<CorrectedMergedRead> {
        // Pull out the ID and the sequence from prior to correction, as we'll be recycling these.
        let id = self.id;
        let seq = self.consensus_seq;

        // Run correction on the quality scores, for which we'll use a handy parallel iterator from rayon
        let corrected_quals = izip!(
                self.fwd_source_seq,
                self.rev_source_seq,
                self.fwd_source_qual,
                self.rev_source_qual,
            )
            .par_bridge() // rust is seriously magic sometimes
            .map(|(fwd_base, rev_base, fwd_qual, rev_qual)| {
                // fill the necessary information for this vertical slice of the overlap into a
                // structure for score correction
                let base_overlap = BaseOverlap {
                    fwd_base,
                    rev_base: &rev_base,
                    fwd_qual,
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
    use crate::merge::UncorrectedMergedRead;

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
        let uncorrected = UncorrectedMergedRead {
            id: "read1".to_string(),
            consensus_seq: b"ACGT".to_vec(),
            consensus_qual: b"IIII".to_vec(),
            fwd_source_seq: b"ACGT",
            fwd_source_qual: b"IIII",
            rev_source_seq: b"ACGT".to_vec(),
            rev_source_qual: b"IIII".to_vec(),
        };

        let corrected = uncorrected.correct().unwrap();
        assert_eq!(corrected.id(), "read1");
        assert_eq!(corrected.sequence(), b"ACGT");
        assert_eq!(corrected.sequence().len(), corrected.qualities().len());
    }

    #[test]
    #[ignore = "Known issue: correction only emits overlap-length qualities"]
    fn test_corrected_qualities_match_consensus_len_with_overhangs() {
        let uncorrected = UncorrectedMergedRead {
            id: "read1".to_string(),
            consensus_seq: b"TTTTACGT".to_vec(),
            consensus_qual: b"IIIIIIII".to_vec(),
            fwd_source_seq: b"ACGT",
            fwd_source_qual: b"IIII",
            rev_source_seq: b"ACGT".to_vec(),
            rev_source_qual: b"IIII".to_vec(),
        };

        let corrected = uncorrected.correct().unwrap();
        assert_eq!(corrected.sequence().len(), corrected.qualities().len());
    }
}
