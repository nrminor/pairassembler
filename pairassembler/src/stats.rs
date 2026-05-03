use std::fmt;

use serde::Serialize;

/// Expected non-merged outcome for a correctly paired input record group.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum UnmergedReason {
    /// No overlap satisfied the configured overlap search thresholds.
    NoAcceptableOverlap,
    /// An overlap was found but rejected by validation.
    OverlapRejectedByValidation,
}

impl fmt::Display for UnmergedReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoAcceptableOverlap => f.write_str("no_acceptable_overlap"),
            Self::OverlapRejectedByValidation => f.write_str("overlap_rejected_by_validation"),
        }
    }
}

/// Running counters for one `pairasm` invocation.
#[derive(Clone, Debug, Default, Serialize)]
pub struct AssemblyStats {
    /// Complete R1/R2 record groups read from input.
    pub pairs_seen: u64,
    /// Correctly paired record groups attempted by the assembly pipeline.
    pub pairs_processed: u64,
    /// Pairs merged into a consensus read.
    pub pairs_merged: u64,
    /// Correctly paired records that were not merged for expected biological/configuration reasons.
    pub pairs_unmerged: u64,
    /// Unmerged pairs with no overlap satisfying search thresholds.
    pub pairs_unmerged_no_overlap: u64,
    /// Unmerged pairs with overlap evidence rejected by validation.
    pub pairs_unmerged_validation_rejected: u64,
    /// Unmerged pairs written to the optional unmerged output.
    pub unmerged_pairs_written: u64,
    /// Unmerged pairs not written because no unmerged output was requested.
    pub unmerged_pairs_not_written: u64,
    /// Input record groups skipped because R1/R2 keys did not agree.
    pub mate_id_mismatches: u64,
    /// Total input bases observed across both mates.
    pub bases_in: u64,
    /// Total consensus bases emitted in merged reads.
    pub bases_merged: u64,
    /// Whether overlap-based quality correction was enabled for this run.
    pub quality_correction_enabled: bool,
}

impl AssemblyStats {
    /// Create counters with run-level configuration recorded for summaries.
    #[must_use]
    pub const fn new(quality_correction_enabled: bool) -> Self {
        Self {
            quality_correction_enabled,
            pairs_seen: 0,
            pairs_processed: 0,
            pairs_merged: 0,
            pairs_unmerged: 0,
            pairs_unmerged_no_overlap: 0,
            pairs_unmerged_validation_rejected: 0,
            unmerged_pairs_written: 0,
            unmerged_pairs_not_written: 0,
            mate_id_mismatches: 0,
            bases_in: 0,
            bases_merged: 0,
        }
    }

    /// Record one complete R1/R2 group read from input.
    pub fn record_pair_seen(&mut self, r1_bases: usize, r2_bases: usize) {
        self.pairs_seen = self.pairs_seen.saturating_add(1);
        self.bases_in = self
            .bases_in
            .saturating_add(usize_to_u64(r1_bases).saturating_add(usize_to_u64(r2_bases)));
    }

    /// Record one input contract violation where R1/R2 keys do not agree.
    pub fn record_mate_id_mismatch(&mut self) {
        self.mate_id_mismatches = self.mate_id_mismatches.saturating_add(1);
    }

    /// Record one successfully merged pair.
    pub fn record_merged(&mut self, merged_bases: usize) {
        self.pairs_processed = self.pairs_processed.saturating_add(1);
        self.pairs_merged = self.pairs_merged.saturating_add(1);
        self.bases_merged = self.bases_merged.saturating_add(usize_to_u64(merged_bases));
    }

    /// Record one correctly paired but unmerged pair.
    pub fn record_unmerged(&mut self, reason: UnmergedReason, was_written: bool) {
        self.pairs_processed = self.pairs_processed.saturating_add(1);
        self.pairs_unmerged = self.pairs_unmerged.saturating_add(1);

        match reason {
            UnmergedReason::NoAcceptableOverlap => {
                self.pairs_unmerged_no_overlap = self.pairs_unmerged_no_overlap.saturating_add(1);
            },
            UnmergedReason::OverlapRejectedByValidation => {
                self.pairs_unmerged_validation_rejected =
                    self.pairs_unmerged_validation_rejected.saturating_add(1);
            },
        }

        if was_written {
            self.unmerged_pairs_written = self.unmerged_pairs_written.saturating_add(1);
        } else {
            self.unmerged_pairs_not_written = self.unmerged_pairs_not_written.saturating_add(1);
        }
    }

    /// Return the most common unmerged reason for progress reporting.
    #[must_use]
    pub fn top_unmerged_reason(&self) -> Option<(UnmergedReason, u64)> {
        let no_overlap = self.pairs_unmerged_no_overlap;
        let validation = self.pairs_unmerged_validation_rejected;

        match (no_overlap, validation) {
            (0, 0) => None,
            (left, right) if left >= right => Some((UnmergedReason::NoAcceptableOverlap, left)),
            (_, right) => Some((UnmergedReason::OverlapRejectedByValidation, right)),
        }
    }
}

fn usize_to_u64(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}

#[cfg(test)]
mod tests {
    use super::{AssemblyStats, UnmergedReason};

    #[test]
    fn stats_account_for_unmerged_reasons_neutrally() {
        let mut stats = AssemblyStats::new(true);

        stats.record_pair_seen(4, 4);
        stats.record_unmerged(UnmergedReason::NoAcceptableOverlap, true);
        stats.record_pair_seen(4, 4);
        stats.record_unmerged(UnmergedReason::OverlapRejectedByValidation, false);

        assert_eq!(stats.pairs_seen, 2);
        assert_eq!(stats.pairs_processed, 2);
        assert_eq!(stats.pairs_unmerged, 2);
        assert_eq!(stats.unmerged_pairs_written, 1);
        assert_eq!(stats.unmerged_pairs_not_written, 1);
    }

    #[test]
    fn stats_do_not_count_mate_mismatch_as_unmerged() {
        let mut stats = AssemblyStats::new(false);

        stats.record_pair_seen(4, 4);
        stats.record_mate_id_mismatch();

        assert_eq!(stats.pairs_seen, 1);
        assert_eq!(stats.mate_id_mismatches, 1);
        assert_eq!(stats.pairs_unmerged, 0);
        assert_eq!(stats.pairs_processed, 0);
    }
}
