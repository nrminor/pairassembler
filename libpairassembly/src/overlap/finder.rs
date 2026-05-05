use crate::Result;
#[cfg(test)]
use crate::read::ReadPair;
use wide::{CmpEq, u8x32};

#[cfg(test)]
use super::AssemblyScratch;
use super::{
    HasOrientedPairSlices, OrientedPairSlices, OverlapBounds, OverlapParams, OverlapSpan,
    PairOverlap,
};

const SIMD_LANES: usize = 32;

/// Finds pair overlaps using configured no-gap overlap search heuristics.
pub(crate) struct OverlapFinder<'params> {
    params: &'params OverlapParams,
}

impl<'params> OverlapFinder<'params> {
    pub(crate) fn new(params: &'params OverlapParams) -> Self {
        Self { params }
    }

    #[cfg(test)]
    pub(crate) fn find<'pair, 'scratch>(
        &self,
        pair: ReadPair<'pair>,
        scratch: &'scratch mut AssemblyScratch,
    ) -> Result<Option<PairOverlap<'pair, 'scratch>>> {
        let slices = pair.to_oriented_slices(scratch);
        self.find_in_slices(slices)
    }

    pub(crate) fn find_in_slices<'pair, 'scratch>(
        &self,
        slices: OrientedPairSlices<'pair, 'scratch>,
    ) -> Result<Option<PairOverlap<'pair, 'scratch>>> {
        let Some(overlap_span) = self.scan_for_overlap_span_both(&slices)? else {
            return Ok(None);
        };

        PairOverlap::from_span(slices, overlap_span).map(Some)
    }

    /// Scan both directional overlap layouts and reconcile them via tie policy.
    ///
    /// # Errors
    ///
    /// Returns tie-policy or overlap-span construction errors from downstream
    /// directional scanners.
    fn scan_for_overlap_span_both(
        &self,
        slices: &OrientedPairSlices<'_, '_>,
    ) -> Result<Option<OverlapSpan>> {
        let (read1, read2) = slices.sequences();
        let overlap_from_left = self.scan_from_start(read1, read2)?;
        let overlap_from_right = self.scan_from_end(read1, read2)?;

        self.resolve_overlap_tie(overlap_from_left, overlap_from_right)
    }

    /// Scan candidate overlaps by sliding the forward read start across read1.
    ///
    /// In this mode, `r2_start` is fixed at 0 while `r1_start` increases.
    fn scan_from_start(&self, read1: &[u8], read2: &[u8]) -> Result<Option<OverlapSpan>> {
        let upper = read1.len().saturating_sub(self.min_overlap());

        for offset in 0..=upper {
            let overlap_len = (read1.len() - offset).min(read2.len());
            if overlap_len < self.min_overlap() {
                break;
            }

            let candidate = Candidate {
                overlap_len,
                r1_start: offset,
                r2_start: 0,
            };

            if let Some(hit) = candidate.evaluate(read1, read2, self)? {
                return Ok(Some(hit));
            }
        }

        Ok(None)
    }

    /// Scan candidate overlaps by sliding the reverse-complemented read start across read2.
    ///
    /// In this mode, `r1_start` is fixed at 0 while `r2_start` increases.
    fn scan_from_end(&self, read1: &[u8], read2: &[u8]) -> Result<Option<OverlapSpan>> {
        // FromEnd semantics (aligned with the oracle/fastp-style no-gap loop):
        //
        //   r1 window: read1[0 .. overlap_len)
        //   r2 window: read2[offset .. offset + overlap_len)
        //
        // In this mode, only the read2 window start shifts with `offset`; read1 start remains 0.
        let upper = read2.len().saturating_sub(self.min_overlap());

        for offset in 0..=upper {
            let overlap_len = read1.len().min(read2.len() - offset);
            if overlap_len < self.min_overlap() {
                break;
            }

            let candidate = Candidate {
                overlap_len,
                r1_start: 0,
                r2_start: offset,
            };

            if let Some(hit) = candidate.evaluate(read1, read2, self)? {
                return Ok(Some(hit));
            }
        }

        Ok(None)
    }

    fn overlap_params(&self) -> &OverlapParams {
        self.params
    }

    fn min_overlap(&self) -> usize {
        self.overlap_params().min_overlap()
    }

    fn min_comparisons(&self) -> usize {
        self.overlap_params().min_comparisons()
    }

    fn overlap_diff_max_for(&self, overlap_len: usize) -> usize {
        self.overlap_params().allowed_differences_for(overlap_len)
    }

    fn resolve_overlap_tie(
        &self,
        from_start_hit: Option<OverlapSpan>,
        from_end_hit: Option<OverlapSpan>,
    ) -> Result<Option<OverlapSpan>> {
        self.overlap_params()
            .tie_policy()
            .resolve(from_start_hit, from_end_hit)
    }
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    overlap_len: usize,
    r1_start: usize,
    r2_start: usize,
}

impl Candidate {
    /// Return borrowed windows over the read pair corresponding to this candidate overlap.
    ///
    /// This keeps candidate-to-sequence mapping explicit while avoiding temporary allocations.
    #[inline]
    fn overlap_windows<'a>(self, read1: &'a [u8], read2: &'a [u8]) -> (&'a [u8], &'a [u8]) {
        let left = &read1[self.r1_start..self.r1_start + self.overlap_len];
        let right = &read2[self.r2_start..self.r2_start + self.overlap_len];
        (left, right)
    }

    #[inline]
    fn evaluate(
        self,
        read1: &[u8],
        read2: &[u8],
        finder: &OverlapFinder<'_>,
    ) -> Result<Option<OverlapSpan>> {
        if self.overlap_len < finder.min_overlap() {
            return Ok(None);
        }

        let (left, right) = self.overlap_windows(read1, read2);
        let overlap_diff_max = finder.overlap_diff_max_for(self.overlap_len);

        let scan = count_mismatches_bounded_simd(left, right, self.overlap_len, overlap_diff_max);
        if scan.exceeded_limit() {
            return Ok(None);
        }
        let diff = scan.mismatches();

        if diff > overlap_diff_max || self.overlap_len < finder.min_comparisons() {
            return Ok(None);
        }

        let bounds = OverlapBounds::new(self.overlap_len, self.r1_start, self.r2_start);
        OverlapSpan::new(bounds, diff, read1.len(), read2.len()).map(Some)
    }
}

#[derive(Debug, Clone, Copy)]
struct MismatchScan {
    mismatch_count: usize,
    exceeded_limit: bool,
}

impl MismatchScan {
    #[inline]
    fn within_limit(mismatch_count: usize) -> Self {
        Self {
            mismatch_count,
            exceeded_limit: false,
        }
    }

    #[inline]
    fn exceeded(mismatch_count: usize) -> Self {
        Self {
            mismatch_count,
            exceeded_limit: true,
        }
    }

    #[inline]
    fn mismatches(self) -> usize {
        self.mismatch_count
    }

    #[inline]
    fn exceeded_limit(self) -> bool {
        self.exceeded_limit
    }
}

#[inline]
fn count_mismatches_bounded_simd(
    left: &[u8],
    right: &[u8],
    compare_len: usize,
    mismatch_limit: usize,
) -> MismatchScan {
    let compare_len = left.len().min(right.len()).min(compare_len);

    let mut mismatches = 0usize;

    let (left_chunks, left_tail) = left[..compare_len].as_chunks::<SIMD_LANES>();
    let (right_chunks, right_tail) = right[..compare_len].as_chunks::<SIMD_LANES>();

    debug_assert_eq!(left_tail.len(), right_tail.len());

    for idx in 0..left_chunks.len() {
        let left_vec = u8x32::from(left_chunks[idx]);
        let right_vec = u8x32::from(right_chunks[idx]);
        let equal_mask = left_vec.simd_eq(right_vec).to_bitmask();
        mismatches += SIMD_LANES - equal_mask.count_ones() as usize;

        if mismatches > mismatch_limit {
            return MismatchScan::exceeded(mismatches);
        }
    }

    for idx in 0..left_tail.len() {
        mismatches += usize::from(left_tail[idx] != right_tail[idx]);
        if mismatches > mismatch_limit {
            return MismatchScan::exceeded(mismatches);
        }
    }

    MismatchScan::within_limit(mismatches)
}

#[cfg(test)]
mod tests {
    use crate::{
        Error,
        errors::OverlapError,
        read::{ReadPair, SequenceRead},
    };
    use proptest::{collection::vec, prelude::*};

    use super::*;
    use crate::overlap::{AssemblyScratch, TiePolicy, slices::reverse_complement_bytes};

    /// Test-only reference implementation that mirrors fastp's simple two-loop no-gap overlap
    /// search. This is intentionally explicit and non-clever so it can serve as a behavioral
    /// oracle for scanner behavior.
    fn oracle_scan_no_gap_from_start(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Option<OverlapSpan> {
        let r1 = mates.fwd_sequence_bytes();
        let r2_rc = reverse_complement_bytes(mates.rev_sequence_bytes());
        let r2 = r2_rc.as_ref();

        let len1 = r1.len();
        let len2 = r2.len();
        let upper = len1.saturating_sub(params.min_overlap());

        for offset in 0..=upper {
            let overlap_len = (len1 - offset).min(len2);
            if overlap_len < params.min_overlap() {
                break;
            }

            let overlap_diff_max = params.allowed_differences_for(overlap_len);

            let mut diff = 0;
            let mut compared = 0;

            for i in 0..overlap_len {
                compared = i + 1;

                if r1[offset + i] != r2[i] {
                    diff += 1;
                    if diff > overlap_diff_max && compared < params.min_comparisons() {
                        break;
                    }
                }
            }

            if diff <= overlap_diff_max && compared >= params.min_comparisons() {
                let bounds = OverlapBounds::new(overlap_len, offset, 0);
                return OverlapSpan::new(bounds, diff, len1, len2).ok();
            }
        }

        None
    }

    fn oracle_scan_no_gap_from_end(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Option<OverlapSpan> {
        let r1 = mates.fwd_sequence_bytes();
        let r2_rc = reverse_complement_bytes(mates.rev_sequence_bytes());
        let r2 = r2_rc.as_ref();

        let len1 = r1.len();
        let len2 = r2.len();
        let upper = len2.saturating_sub(params.min_overlap());

        for k in 0..=upper {
            let overlap_len = len1.min(len2 - k);
            if overlap_len < params.min_overlap() {
                break;
            }

            let overlap_diff_max = params.allowed_differences_for(overlap_len);

            let mut diff = 0;
            let mut compared = 0;

            for i in 0..overlap_len {
                compared = i + 1;

                if r1[i] != r2[k + i] {
                    diff += 1;
                    if diff > overlap_diff_max && compared < params.min_comparisons() {
                        break;
                    }
                }
            }

            if diff <= overlap_diff_max && compared >= params.min_comparisons() {
                let bounds = OverlapBounds::new(overlap_len, 0, k);
                return OverlapSpan::new(bounds, diff, len1, len2).ok();
            }
        }

        None
    }

    fn oracle_scan_no_gap(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Result<Option<OverlapSpan>> {
        let from_start = oracle_scan_no_gap_from_start(mates, params);
        let from_end = oracle_scan_no_gap_from_end(mates, params);
        params.tie_policy().resolve(from_start, from_end)
    }

    fn scan_bounds(mates: &ReadPair<'_>, params: &OverlapParams) -> Result<Option<OverlapSpan>> {
        let mut scratch = AssemblyScratch::default();
        let slices = mates.to_oriented_slices(&mut scratch);
        OverlapFinder::new(params).scan_for_overlap_span_both(&slices)
    }

    fn scan_bounds_from_start(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Result<Option<OverlapSpan>> {
        let mut scratch = AssemblyScratch::default();
        let slices = mates.to_oriented_slices(&mut scratch);
        let finder = OverlapFinder::new(params);
        let (read1, read2) = slices.sequences();
        finder.scan_from_start(read1, read2)
    }

    fn scan_bounds_from_end(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Result<Option<OverlapSpan>> {
        let mut scratch = AssemblyScratch::default();
        let slices = mates.to_oriented_slices(&mut scratch);
        let finder = OverlapFinder::new(params);
        let (read1, read2) = slices.sequences();
        finder.scan_from_end(read1, read2)
    }

    #[test]
    fn test_scan_matches_oracle_from_start_simple_case() {
        let r1 = SequenceRead::new("read1", "TTTTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let params = OverlapParams::default().with_settings(2, 4, 0.2, 4);

        let expected = oracle_scan_no_gap_from_start(&mates, &params);
        let observed = scan_bounds_from_start(&mates, &params)
            .expect("from-start scanner should not error on simple oracle fixture");

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_scan_matches_oracle_from_end_simple_case() {
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGTAAAA", "IIIIIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let params = OverlapParams::default().with_settings(2, 4, 0.2, 4);

        let expected = oracle_scan_no_gap_from_end(&mates, &params);
        let observed = scan_bounds_from_end(&mates, &params)
            .expect("from-end scanner should not error on simple oracle fixture");

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_overlap_from_start_correct_bounds() {
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGT", "IIII");

        let params = OverlapParams::default().with_settings(2, 3, 0.2, 3);

        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let overlap = scan_bounds_from_start(&mates, &params)
            .expect("from-start scanner should not error when checking canonical bounds");

        assert!(overlap.is_some());
        let span = overlap.expect("expected overlap in canonical from-start bounds fixture");
        assert_eq!(span.bounds().fwd_start_offset(), 0);
        assert_eq!(span.bounds().fwd_end_offset(), 3);
        assert_eq!(span.bounds().rev_start_offset(), 0);
        assert_eq!(span.bounds().rev_end_offset(), 3);
    }

    #[test]
    fn test_overlap_from_end_correct_bounds() {
        let r1 = SequenceRead::new("read1", "TTTTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTAAAA", "IIIIIIII");

        let params = OverlapParams::default().with_settings(1, 4, 0.1, 4);

        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let overlap = scan_bounds_from_end(&mates, &params)
            .expect("from-end scanner should not error when checking canonical bounds");

        assert!(overlap.is_some());
        let span = overlap.expect("expected overlap in canonical from-end bounds fixture");
        assert_eq!(span.bounds().fwd_start_offset(), 0);
        assert_eq!(span.bounds().fwd_end_offset(), 7);
        assert_eq!(span.bounds().rev_start_offset(), 0);
        assert_eq!(span.bounds().rev_end_offset(), 7);
    }

    #[test]
    fn test_no_overlap_detected() {
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "TTTT", "IIII");

        let params = OverlapParams::default().with_settings(1, 4, 0.1, 4);

        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let overlap = scan_bounds_from_start(&mates, &params)
            .expect("from-start scanner should not error in no-overlap fixture");

        assert!(overlap.is_none());
    }

    #[derive(Clone)]
    struct OverlapFixture {
        name: &'static str,
        r1: &'static str,
        r2: &'static str,
        params: OverlapParams,
        expect_overlap: bool,
    }

    #[test]
    fn test_overlap_edge_case_fixtures_match_oracle() {
        let fixtures = vec![
            OverlapFixture {
                name: "exact_min_overlap_at_boundary_is_scanned",
                r1: "GGGGACGT",
                r2: "ACGT",
                params: OverlapParams::default().with_settings(0, 4, 0.0, 4),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "below_min_overlap_rejected",
                r1: "GGGGACGT",
                r2: "ACGTA",
                params: OverlapParams::default().with_settings(0, 5, 0.0, 4),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "diff_equal_threshold_accepted",
                r1: "ACGTACGT",
                r2: "TCGTACGT",
                params: OverlapParams::default().with_settings(1, 7, 0.2, 8),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "diff_above_threshold_rejected",
                r1: "ACGTACGT",
                r2: "TCGTACGT",
                params: OverlapParams::default().with_settings(0, 7, 0.0, 8),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "short_perfect_but_below_min_comparisons_rejected",
                r1: "ACGT",
                r2: "ACGT",
                params: OverlapParams::default().with_settings(0, 4, 0.0, 5),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "asymmetric_lengths_from_end_detected",
                r1: "ACGTACGT",
                r2: "ACGTACGTAAAA",
                params: OverlapParams::default().with_settings(2, 4, 0.2, 4),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "low_complexity_no_overlap",
                r1: "AAAAAAAAAAAA",
                r2: "CCCCCCCCCCCC",
                params: OverlapParams::default().with_settings(0, 6, 0.0, 6),
                expect_overlap: false,
            },
        ];

        for fixture in fixtures {
            let q1 = "I".repeat(fixture.r1.len());
            let q2 = "I".repeat(fixture.r2.len());
            let r1 = SequenceRead::new("read1", fixture.r1, &q1);
            let r2 = SequenceRead::new("read1", fixture.r2, &q2);
            let mates =
                ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

            let expected = oracle_scan_no_gap(&mates, &fixture.params).unwrap_or_else(|err| {
                panic!(
                    "fixture '{}' oracle errored unexpectedly: {err}",
                    fixture.name
                )
            });
            let observed = scan_bounds(&mates, &fixture.params).unwrap_or_else(|err| {
                panic!("fixture '{}' errored unexpectedly: {err}", fixture.name)
            });
            assert_eq!(
                observed, expected,
                "fixture '{}' diverged from oracle",
                fixture.name
            );
            assert_eq!(
                observed.is_some(),
                fixture.expect_overlap,
                "fixture '{}' overlap expectation mismatch",
                fixture.name
            );
        }
    }

    #[test]
    fn test_scan_matches_oracle_from_start_simple_case_canonical() {
        let r1 = SequenceRead::new("read1", "TTTTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let params = OverlapParams::default().with_settings(2, 4, 0.2, 4);

        let oracle = oracle_scan_no_gap_from_start(&mates, &params);
        let got = scan_bounds_from_start(&mates, &params)
            .expect("from-start scanner should not error for canonical oracle test");

        assert_eq!(got, oracle);
    }

    #[test]
    fn test_scan_matches_oracle_from_end_simple_case_canonical() {
        let r1 = SequenceRead::new("read1", "TTTTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let params = OverlapParams::default().with_settings(2, 4, 0.2, 4);

        let oracle = oracle_scan_no_gap_from_end(&mates, &params);
        let got = scan_bounds_from_end(&mates, &params)
            .expect("from-end scanner should not error for canonical oracle test");

        assert_eq!(got, oracle);
    }

    #[test]
    fn test_scan_matches_oracle_regression_from_end_case_canonical() {
        let r1 = SequenceRead::new("read1", "AGTGAAGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "TCACAAAAA", "IIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let params = OverlapParams::new(1, 4, 0.296_884_83, 1);

        let oracle = oracle_scan_no_gap_from_end(&mates, &params);
        let got = scan_bounds_from_end(&mates, &params)
            .expect("from-end scanner should not error for regression oracle test");

        assert_eq!(got, oracle);
    }

    #[test]
    fn test_from_end_shifts_only_r2_window_start() {
        let r1 = SequenceRead::new("read1", "AGTGAAGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "TCACAAAAA", "IIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let params = OverlapParams::new(1, 4, 0.296_884_83, 1);

        let got = scan_bounds_from_end(&mates, &params)
            .expect("from-end scanner should not error for window-shift regression")
            .expect("expected overlap in from-end window-shift regression");

        assert_eq!(got.bounds().fwd_start_offset(), 0);
        assert_eq!(got.bounds().rev_start_offset(), 4);
        assert_eq!(got.overlap_len(), 5);
        assert_eq!(got.diff(), 1);
    }

    #[test]
    fn test_overlap_rejects_direction_tie() {
        let r1 = SequenceRead::new("read-tie", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read-tie", "ACGTACGT", "IIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let params = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);

        let mut scratch = AssemblyScratch::default();
        let got = mates.overlap(&params, &mut scratch);
        assert!(matches!(
            got,
            Err(Error::Overlap(OverlapError::OverlapTie { .. }))
        ));
    }

    fn dna_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        vec(
            prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
            min_len..=max_len,
        )
        .prop_map(|chars| chars.into_iter().collect())
    }

    proptest! {
        #[test]
        fn proptest_scan_matches_oracle_no_gap(
            r1 in dna_string_strategy(8, 80),
            r2 in dna_string_strategy(8, 80),
            overlap_diff_max in 0usize..=4,
            diff_percent_max in 0.0f32..=0.35f32,
            min_overlap_raw in 4usize..=32,
            min_comparisons_raw in 1usize..=64,
        ) {
            let max_possible_overlap = r1.len().min(r2.len());
            let min_overlap = min_overlap_raw.min(max_possible_overlap).max(1);
            let min_comparisons = min_comparisons_raw.min(max_possible_overlap).max(1);

            let params = OverlapParams::new(
                overlap_diff_max,
                min_overlap,
                diff_percent_max,
                min_comparisons,
            );

            let q1 = "I".repeat(r1.len());
            let q2 = "I".repeat(r2.len());
            let mates = ReadPair::from(
                SequenceRead::new("read1", &r1, &q1),
                SequenceRead::new("read1", &r2, &q2),
            )
            .expect("proptest fixture reads should share the same id");

            let observed = scan_bounds(&mates, &params);
            let expected = oracle_scan_no_gap(&mates, &params);
            prop_assert!(observed.is_ok(), "scanner returned unexpected error: {observed:?}");
            prop_assert!(expected.is_ok(), "oracle returned unexpected error: {expected:?}");

            let observed = observed
                .expect("scanner should not error for proptest-generated input");
            let expected = expected
                .expect("oracle should not error for proptest-generated input");
            prop_assert_eq!(&observed, &expected, "scanner mismatch");

            if let Some(hit) = observed {
                let bounds = hit.bounds();
                prop_assert!(bounds.fwd_start_offset() <= bounds.fwd_end_offset());
                prop_assert!(bounds.rev_start_offset() <= bounds.rev_end_offset());
                prop_assert!(bounds.fwd_end_offset() < r1.len());
                prop_assert!(bounds.rev_end_offset() < r2.len());
                prop_assert_eq!(bounds.fwd_end_offset() - bounds.fwd_start_offset() + 1, hit.overlap_len());
                prop_assert_eq!(bounds.rev_end_offset() - bounds.rev_start_offset() + 1, hit.overlap_len());
            }

        }
    }
}
