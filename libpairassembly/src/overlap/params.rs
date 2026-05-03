use crate::{Result, errors::OverlapError};

use super::OverlapSpan;

/// Parameters controlling no-gap overlap search thresholds.
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// Minimum number of base comparisons required before an overlap can be accepted.
    min_comparisons: usize,
    tie_policy: TiePolicy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Policy for handling equal-quality overlaps found from both search directions.
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
    /// The candidate with lower mismatch rate wins (`diff / overlap_len`).
    /// Exact-rate ties are handled according to the selected policy.
    ///
    /// # Errors
    ///
    /// Returns `OverlapTie` when both candidates have equal mismatch rate and
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
                let left_key = left.diff() * right.overlap_len();
                let right_key = right.diff() * left.overlap_len();

                if left_key < right_key {
                    return Ok(Some(left));
                }
                if right_key < left_key {
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

impl Default for OverlapParams {
    fn default() -> Self {
        OverlapParams {
            overlap_diff_max: 2,
            min_overlap: 30,
            diff_percent_max: 0.2,
            min_comparisons: 50,
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
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn allowed_differences_for(&self, overlap_len: usize) -> usize {
        self.overlap_diff_max()
            .min((overlap_len as f32 * self.diff_percent_max()) as usize)
    }
}
