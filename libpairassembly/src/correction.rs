#[derive(Debug, Clone)]
pub struct BaseOverlap<'overlap> {
    fwd_base: &'overlap u8,
    rev_base: &'overlap u8,
    fwd_qual: &'overlap i32,
    rev_qual: &'overlap i32,
}

impl<'overlap> BaseOverlap<'overlap> {
    pub fn new(
        fwd_base: &'overlap u8,
        rev_base: &'overlap u8,
        fwd_qual: &'overlap i32,
        rev_qual: &'overlap i32,
    ) -> Self {
        Self {
            fwd_base,
            rev_base,
            fwd_qual,
            rev_qual,
        }
    }

    pub fn compute_overlap_score(&self) -> (&'overlap u8, u8) {
        let fwd_error = 10_f64.powf(-(self.fwd_qual.to_owned() / 10) as f64);
        let rev_error = 10_f64.powf(-(self.rev_qual.to_owned() / 10) as f64);

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

impl<'err_prob> MatchStatus<'err_prob> {
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

fn mismatch_error_probability(fwd_error: &f64, rev_error: &f64) -> f64 {
    ((fwd_error * rev_error) / 3.0)
        / ((1.0 - fwd_error) * (1.0 - rev_error) + 4.0 * (fwd_error * rev_error) / 3.0)
}

fn match_error_probability(fwd_error: &f64, rev_error: &f64) -> f64 {
    (fwd_error * (1.0 - rev_error / 3.0))
        / (fwd_error + rev_error - 4.0 * (fwd_error * rev_error) / 3.0)
}
