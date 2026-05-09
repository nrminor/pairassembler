use crate::{Result, errors::OverlapError};

use super::OverlapSpan;

/// Parameters controlling no-gap overlap search thresholds.
///
/// The overlap finder searches for ungapped suffix/prefix overlaps between mates. An overlap is
/// accepted only when it satisfies both absolute and percent mismatch limits, and has enough
/// base-to-base comparisons for the configured thresholds to be meaningful.
///
/// ```rust
/// use libpairassembly::{OverlapParams, TiePolicy};
///
/// let params = OverlapParams::default()
///     .with_min_overlap(40)
///     .with_tie_policy(TiePolicy::PreferFromStart);
///
/// assert_eq!(params.min_overlap(), 40);
/// assert_eq!(params.allowed_differences_for(50), 5);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// Minimum number of base comparisons required before an overlap can be accepted.
    min_comparisons: usize,
    tie_policy: TiePolicy,
}

/// Policy for handling equal-quality overlaps found from both search directions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TiePolicy {
    /// Return an `OverlapTie` error.
    Reject,
    /// Keep the overlap found from the start-oriented search.
    PreferFromStart,
    /// Keep the overlap found from the end-oriented search.
    PreferFromEnd,
}

impl TiePolicy {
    /// Resolve two directional overlap candidates into a single winner.
    ///
    /// The candidate with stronger overlap evidence wins. The score rewards
    /// overlap length and penalizes mismatches so that short exact overlaps do
    /// not beat much longer near-exact overlaps merely by having a lower
    /// mismatch rate. Exact evidence-score ties are handled according to the
    /// selected policy.
    ///
    /// # Errors
    ///
    /// Returns `OverlapTie` when both candidates have equal evidence score and
    /// the policy is [`TiePolicy::Reject`].
    pub(crate) fn resolve(
        self,
        from_start_hit: Option<OverlapSpan>,
        from_end_hit: Option<OverlapSpan>,
    ) -> Result<Option<OverlapSpan>> {
        match (from_start_hit, from_end_hit) {
            (None, None) => Ok(None),
            (Some(left), None) => Ok(Some(left)),
            (None, Some(right)) => Ok(Some(right)),
            (Some(left), Some(right)) => {
                let left_key = overlap_evidence_score(left);
                let right_key = overlap_evidence_score(right);

                if left_key > right_key {
                    return Ok(Some(left));
                }
                if right_key > left_key {
                    return Ok(Some(right));
                }

                match self {
                    TiePolicy::Reject => Err(OverlapError::OverlapTie {
                        diff: left.diff(),
                        overlap_len: left.overlap_len(),
                    }
                    .into()),
                    TiePolicy::PreferFromStart => Ok(Some(left)),
                    TiePolicy::PreferFromEnd => Ok(Some(right)),
                }
            },
        }
    }
}

fn overlap_evidence_score(overlap: OverlapSpan) -> isize {
    // Count a mismatch as one lost matching base plus an additional four-base
    // penalty. This keeps exact short overlaps competitive with modestly longer
    // noisy overlaps, but prevents tiny exact overlaps from beating much longer
    // near-exact evidence.
    const MISMATCH_PENALTY: isize = 5;

    let overlap_len = isize::try_from(overlap.overlap_len()).unwrap_or(isize::MAX);
    let diff = isize::try_from(overlap.diff()).unwrap_or(isize::MAX);

    overlap_len.saturating_sub(diff.saturating_mul(MISMATCH_PENALTY))
}

impl Default for OverlapParams {
    fn default() -> Self {
        OverlapParams {
            overlap_diff_max: 5,
            min_overlap: 30,
            diff_percent_max: 0.2,
            min_comparisons: 30,
            tie_policy: TiePolicy::PreferFromStart,
        }
    }
}

impl OverlapParams {
    #[must_use]
    pub fn new(
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
    ) -> Self {
        Self {
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            tie_policy: TiePolicy::PreferFromStart,
        }
    }

    #[must_use]
    pub fn with_settings(
        mut self,
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
    ) -> Self {
        self.overlap_diff_max = overlap_diff_max;
        self.min_overlap = min_overlap;
        self.diff_percent_max = diff_percent_max;
        self.min_comparisons = min_comparisons;
        self
    }

    #[must_use]
    pub fn with_overlap_diff_max(mut self, val: usize) -> Self {
        self.overlap_diff_max = val;
        self
    }

    #[must_use]
    pub fn with_min_overlap(mut self, val: usize) -> Self {
        self.min_overlap = val;
        self
    }

    #[must_use]
    pub fn with_diff_percent_max(mut self, val: f32) -> Self {
        self.diff_percent_max = val;
        self
    }

    #[must_use]
    pub fn with_min_comparisons(mut self, val: usize) -> Self {
        self.min_comparisons = val;
        self
    }

    #[must_use]
    pub fn with_tie_policy(mut self, tie_policy: TiePolicy) -> Self {
        self.tie_policy = tie_policy;
        self
    }

    #[must_use]
    pub fn overlap_diff_max(&self) -> usize {
        self.overlap_diff_max
    }

    #[must_use]
    pub fn min_overlap(&self) -> usize {
        self.min_overlap
    }

    #[must_use]
    pub fn diff_percent_max(&self) -> f32 {
        self.diff_percent_max
    }

    #[must_use]
    pub fn min_comparisons(&self) -> usize {
        self.min_comparisons
    }

    #[must_use]
    pub fn tie_policy(&self) -> TiePolicy {
        self.tie_policy
    }

    #[must_use]
    pub fn allowed_differences_for(&self, overlap_len: usize) -> usize {
        let cap = self.overlap_diff_max();
        let scaled_overlap = usize_to_f64(overlap_len) * f64::from(self.diff_percent_max());

        if scaled_overlap.is_nan() || scaled_overlap <= 0.0 {
            return 0;
        }
        if scaled_overlap >= usize_to_f64(cap) {
            return cap;
        }

        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "scaled_overlap is already finite, positive, and below the usize cap; truncating the floor is the intended threshold calculation"
        )]
        let mut allowed = scaled_overlap.floor() as usize;
        if allowed > cap {
            allowed = cap;
        }

        while allowed > 0 && usize_to_f64(allowed) > scaled_overlap {
            allowed -= 1;
        }

        while allowed < cap && usize_to_f64(allowed + 1) <= scaled_overlap {
            allowed += 1;
        }

        allowed
    }
}

fn usize_to_f64(value: usize) -> f64 {
    const U32_RADIX_USIZE: usize = 4_294_967_296;
    const U32_RADIX_F64: f64 = 4_294_967_296.0;

    let high = value / U32_RADIX_USIZE;
    let low = value % U32_RADIX_USIZE;

    let Ok(high) = u32::try_from(high) else {
        return f64::INFINITY;
    };
    let Ok(low) = u32::try_from(low) else {
        unreachable!("value modulo 2^32 always fits in u32")
    };

    f64::from(high) * U32_RADIX_F64 + f64::from(low)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn allowed_differences_reference(
        overlap_len: usize,
        overlap_diff_max: usize,
        diff_percent_max: f32,
    ) -> usize {
        let scaled_overlap = usize_to_f64(overlap_len) * f64::from(diff_percent_max);

        let mut accepted = 0;
        let mut rejected = overlap_diff_max;
        while accepted < rejected {
            let candidate = accepted + (rejected - accepted).div_ceil(2);
            if usize_to_f64(candidate) <= scaled_overlap {
                accepted = candidate;
            } else {
                rejected = candidate - 1;
            }
        }

        accepted
    }

    proptest! {
        #[test]
        fn allowed_differences_matches_previous_binary_search(
            overlap_len in any::<usize>(),
            overlap_diff_max in any::<usize>(),
            diff_percent_max in -10.0f32..=10.0,
        ) {
            let params = OverlapParams::default()
                .with_overlap_diff_max(overlap_diff_max)
                .with_diff_percent_max(diff_percent_max);

            prop_assert_eq!(
                params.allowed_differences_for(overlap_len),
                allowed_differences_reference(overlap_len, overlap_diff_max, diff_percent_max),
            );
        }
    }
}
