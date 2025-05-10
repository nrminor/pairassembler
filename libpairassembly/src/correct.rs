use std::str;

use crate::{
    SequenceRead,
    merge::{Merge, UncorrectedMergedRead},
};

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
    pub fn correct_quality_scores(self) -> color_eyre::Result<CorrectedMergedRead> {
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
            .map(|(fwd_base, rev_base, fwd_qual, rev_qual)| BaseOverlap {
                fwd_base,
                rev_base,
                fwd_qual,
                rev_qual,
            })
            .map(|base_overlap| {
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

    pub fn compute_corrected_score(&self) -> (&'overlap u8, u8) {
        let fwd_qual = *self.fwd_qual as i32;
        let rev_qual = *self.rev_qual as i32;
        let fwd_error = 10_f64.powf(-(fwd_qual.to_owned() / 10) as f64);
        let rev_error = 10_f64.powf(-(rev_qual.to_owned() / 10) as f64);

        match self.fwd_base == self.rev_base {
            // probability that the two matching self represent an error,
            // given their quality scores
            true => {
                let status = Match {
                    fwd_error: &fwd_error,
                    rev_error: &rev_error,
                };
                let score = status.compute_score();
                (self.fwd_base, score)
            },
            // probability that the two mismatching self represent an error,
            // given their quality scores
            false => {
                let status = Mismatch {
                    fwd_error: &fwd_error,
                    rev_error: &rev_error,
                };
                let score = status.compute_score();
                match self.fwd_qual >= self.rev_qual {
                    true => (self.fwd_base, score),
                    false => (self.rev_base, score),
                }
            },
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
use itertools::izip;
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    slice::ParallelSliceMut,
};

impl MatchStatus<'_> {
    pub fn compute_score(self) -> u8 {
        let posterior = match self {
            Match {
                fwd_error,
                rev_error,
            } => mismatch_error_probability(fwd_error, rev_error),
            Mismatch {
                fwd_error,
                rev_error,
            } => match_error_probability(fwd_error, rev_error),
        };

        // compute the integer quality score
        let score = (posterior.log10() * -10.0).floor();

        if score > 40.0 { 40_u8 } else { score as u8 }
    }
}

#[inline]
fn mismatch_error_probability(fwd_error: &f64, rev_error: &f64) -> f64 {
    ((fwd_error * rev_error) / 3.0)
        / ((1.0 - fwd_error) * (1.0 - rev_error) + 4.0 * (fwd_error * rev_error) / 3.0)
}

#[inline]
fn match_error_probability(fwd_error: &f64, rev_error: &f64) -> f64 {
    (fwd_error * (1.0 - rev_error / 3.0))
        / (fwd_error + rev_error - 4.0 * (fwd_error * rev_error) / 3.0)
}
