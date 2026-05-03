use std::ops::Range;

use crate::{
    ReadPair, Result,
    errors::OverlapError::{
        IndexOutOfBounds, InvalidOverlapLength, OrientedPairSequenceQualityLengthMismatch,
        OverlapTie,
    },
    prelude::utils::decode_fastq_quality_scores,
};
use wide::{CmpEq, u8x32};

const SIMD_LANES: usize = 32;

/// Parameters for the overlap analysis, mostly for skipping read pairs that no need further confirmation
/// that no overlap exists
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// set the minimum amount of base comparisons required to determine if two reads overlap
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
    fn resolve(
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
                    TiePolicy::Reject => Err(OverlapTie {
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
    // TODO: Do some research and tweak these defaults as needed; these are basically just the
    // defaults used in `fastp`, which, while justifiable, are simple heuristics. If anything,
    // because these parameters form the "floor" with respect to overlap quality, it might be
    // worth lowering these heuristics so that validation does more of the work given that it's
    // probabilistic...and statistics > heuristics?
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

impl<'a> ReadPair<'a> {
    /// Discover the best overlap between this read pair.
    ///
    /// Returns `Ok(None)` when no overlap candidate satisfies the configured thresholds.
    ///
    /// # Errors
    ///
    /// Returns an error when overlap candidate reconciliation fails (for example, tie rejection),
    /// or when computed overlap coordinates are inconsistent with read bounds.
    pub fn overlap(&self, params: &OverlapParams) -> Result<Option<PairOverlap<'a>>> {
        OverlapFinder::new(params).find(*self)
    }

    pub(crate) fn to_oriented_slices(self) -> OrientedPairSlices<'a> {
        OrientedPairSlices::from_fastq_ascii_parts(
            self.fwd_id(),
            self.fwd_sequence_bytes(),
            self.fwd_quality_bytes(),
            self.rev_sequence_bytes(),
            self.rev_quality_bytes(),
        )
    }
}

/// Finds pair overlaps using configured no-gap overlap search heuristics.
pub(crate) struct OverlapFinder<'params> {
    params: &'params OverlapParams,
}

impl<'params> OverlapFinder<'params> {
    pub(crate) fn new(params: &'params OverlapParams) -> Self {
        Self { params }
    }

    pub(crate) fn find<'pair>(&self, pair: ReadPair<'pair>) -> Result<Option<PairOverlap<'pair>>> {
        let slices = pair.to_oriented_slices();
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
        slices: &OrientedPairSlices<'_>,
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

impl<'a> OrientedPairSlices<'a> {
    fn from_fastq_ascii_parts(
        id: &'a str,
        fwd_seq: &'a [u8],
        fwd_qual: &[u8],
        rev_raw_seq: &'a [u8],
        rev_raw_qual: &[u8],
    ) -> Self {
        let fwd_qual = decode_fastq_quality_scores(fwd_qual);
        let mut rev_qual_rev = decode_fastq_quality_scores(rev_raw_qual);
        rev_qual_rev.reverse();
        let rev_seq_rc = reverse_complement_bytes(rev_raw_seq);

        Self {
            id,
            fwd_seq,
            fwd_qual,
            rev_seq_rc,
            rev_qual_rev,
        }
    }
}

#[derive(Debug, Clone)]
/// Retained pair slices after converting reads into overlap orientation.
///
/// Forward sequence bytes are borrowed from the original read. Quality scores and reverse-mate
/// orientation are owned fixed-size buffers because FASTQ ASCII qualities must be decoded and the
/// reverse mate must be reverse-complemented before downstream validation, merge, and correction.
pub(crate) struct OrientedPairSlices<'a> {
    pub(crate) id: &'a str,
    pub(crate) fwd_seq: &'a [u8],
    pub(crate) fwd_qual: Box<[u8]>,
    pub(crate) rev_seq_rc: Box<[u8]>,
    pub(crate) rev_qual_rev: Box<[u8]>,
}

/// Exposes score-space paired slices in overlap orientation.
///
/// Implementors must preserve these invariants:
/// - forward sequence is in forward-read orientation;
/// - forward qualities are numeric quality score bytes;
/// - reverse sequence is reverse-complemented into forward-read orientation;
/// - reverse qualities are reversed to match the reverse-complemented sequence;
/// - sequence and quality lengths match within each mate.
pub(crate) trait HasOrientedPairSlices: private::Sealed {
    fn pair_id(&self) -> &str;
    fn forward_sequence(&self) -> &[u8];
    fn forward_quality_score_bytes(&self) -> &[u8];
    fn reverse_sequence_rc(&self) -> &[u8];
    fn reverse_quality_score_bytes_rc(&self) -> &[u8];

    fn forward_len(&self) -> usize {
        self.forward_sequence().len()
    }

    fn reverse_len(&self) -> usize {
        self.reverse_sequence_rc().len()
    }

    fn sequences(&self) -> (&[u8], &[u8]) {
        (self.forward_sequence(), self.reverse_sequence_rc())
    }

    fn quality_score_bytes(&self) -> (&[u8], &[u8]) {
        (
            self.forward_quality_score_bytes(),
            self.reverse_quality_score_bytes_rc(),
        )
    }

    fn validate_shape(&self) -> Result<()> {
        ensure_sequence_quality_lengths(
            "fwd_mate",
            self.forward_sequence().len(),
            self.forward_quality_score_bytes().len(),
        )?;
        ensure_sequence_quality_lengths(
            "rev_mate_rc",
            self.reverse_sequence_rc().len(),
            self.reverse_quality_score_bytes_rc().len(),
        )
    }

    fn validate_overlap_bounds(&self, bounds: OverlapBounds) -> Result<()> {
        self.validate_shape()?;

        if bounds.overlap_len() == 0 {
            return Err(InvalidOverlapLength {
                computed: bounds.overlap_len(),
                read1_len: self.forward_len(),
                read2_len: self.reverse_len(),
                min_required: 1,
            }
            .into());
        }

        let fwd_end = bounds.fwd_end_offset();
        let rev_end = bounds.rev_end_offset();

        let fwd_len = self.forward_len();
        if fwd_end >= fwd_len {
            return Err(IndexOutOfBounds {
                read: "fwd_mate",
                index: fwd_end,
                length: fwd_len,
            }
            .into());
        }

        let rev_len = self.reverse_len();
        if rev_end >= rev_len {
            return Err(IndexOutOfBounds {
                read: "rev_mate",
                index: rev_end,
                length: rev_len,
            }
            .into());
        }

        Ok(())
    }
}

pub(crate) mod private {
    pub(crate) trait Sealed {}
}

impl private::Sealed for OrientedPairSlices<'_> {}

impl HasOrientedPairSlices for OrientedPairSlices<'_> {
    fn pair_id(&self) -> &str {
        self.id
    }

    fn forward_sequence(&self) -> &[u8] {
        self.fwd_seq
    }

    fn forward_quality_score_bytes(&self) -> &[u8] {
        &self.fwd_qual
    }

    fn reverse_sequence_rc(&self) -> &[u8] {
        &self.rev_seq_rc
    }

    fn reverse_quality_score_bytes_rc(&self) -> &[u8] {
        &self.rev_qual_rev
    }
}

fn ensure_sequence_quality_lengths(
    mate: &'static str,
    seq_len: usize,
    qual_len: usize,
) -> Result<()> {
    if seq_len != qual_len {
        return Err(OrientedPairSequenceQualityLengthMismatch {
            mate,
            seq_len,
            qual_len,
        }
        .into());
    }

    Ok(())
}

fn reverse_complement_bytes(seq: &[u8]) -> Box<[u8]> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
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

        // here we set the bounds of the view into the two candidate mates' sequences
        // for this comparison and also compute the max number of allowed differences.
        // Remember: whole read sequences are passed to this function, but its role is
        // only to compare within overlap windows that are static for this function's
        // execution.
        let (left, right) = self.overlap_windows(read1, read2);
        let overlap_diff_max = finder.overlap_diff_max_for(self.overlap_len);

        // Run a single bounded mismatch scan across the whole overlap window.
        // The scan short-circuits as soon as mismatch_limit is exceeded.
        let scan = count_mismatches_bounded_simd(left, right, self.overlap_len, overlap_diff_max);
        if scan.exceeded_limit() {
            return Ok(None);
        }
        let diff = scan.mismatches();

        // by now, we've evaluated the overlap between windows of the two reads set earlier
        // in this function. It's possible that that overlap may have too many differences,
        // or may be too short an overlap. In this case, return an Ok of a None--the
        // search ran into no errors and just came up empty-handed.
        if diff > overlap_diff_max || self.overlap_len < finder.min_comparisons() {
            return Ok(None);
        }

        // otherwise, we have an overlap that's good enough to pass along!
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
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
fn compute_overlap_diff_max(overlap_len: usize, params: &OverlapParams) -> usize {
    params.allowed_differences_for(overlap_len)
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

    // remember that the slices given to this function were trimmed with `Candidate::overlap_windows()`
    // to be the same length, so regardless of where an offset has moved in the loop below,
    // the scalar tails of each read will be the same length
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

/// Canonical overlap span representation used by overlap scanning internals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OverlapSpan {
    bounds: OverlapBounds,
    diff: usize,
}

impl OverlapSpan {
    fn new(bounds: OverlapBounds, diff: usize, r1_len: usize, r2_len: usize) -> Result<Self> {
        let overlap_len = bounds.overlap_len();
        if overlap_len == 0 {
            return Err(InvalidOverlapLength {
                computed: overlap_len,
                read1_len: r1_len,
                read2_len: r2_len,
                min_required: 1,
            }
            .into());
        }

        let r1_end = bounds.forward_range().end;
        let r2_end = bounds.reverse_range().end;

        if r1_end > r1_len {
            return Err(IndexOutOfBounds {
                read: "r1",
                index: r1_end - 1,
                length: r1_len,
            }
            .into());
        }

        if r2_end > r2_len {
            return Err(IndexOutOfBounds {
                read: "r2",
                index: r2_end - 1,
                length: r2_len,
            }
            .into());
        }

        if diff > overlap_len {
            return Err(InvalidOverlapLength {
                computed: diff,
                read1_len: overlap_len,
                read2_len: overlap_len,
                min_required: 0,
            }
            .into());
        }

        Ok(Self { bounds, diff })
    }

    #[inline]
    fn bounds(self) -> OverlapBounds {
        self.bounds
    }

    #[inline]
    fn overlap_len(self) -> usize {
        self.bounds.overlap_len()
    }

    #[inline]
    fn diff(self) -> usize {
        self.diff
    }

    #[inline]
    fn r1_end_inclusive(&self) -> usize {
        self.bounds.fwd_end_offset()
    }

    #[inline]
    fn r2_end_inclusive(&self) -> usize {
        self.bounds.rev_end_offset()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OverlapBounds {
    overlap_len: usize,
    r1_start_offset: usize,
    r2_start_offset: usize,
}

impl OverlapBounds {
    #[inline]
    pub(crate) fn new(overlap_len: usize, r1_start_offset: usize, r2_start_offset: usize) -> Self {
        Self {
            overlap_len,
            r1_start_offset,
            r2_start_offset,
        }
    }

    #[inline]
    pub(crate) fn overlap_len(self) -> usize {
        self.overlap_len
    }

    #[inline]
    pub(crate) fn fwd_start_offset(self) -> usize {
        self.r1_start_offset
    }

    #[inline]
    pub(crate) fn fwd_end_offset(self) -> usize {
        self.r1_start_offset + self.overlap_len - 1
    }

    #[inline]
    pub(crate) fn forward_range(self) -> Range<usize> {
        self.r1_start_offset..self.r1_start_offset + self.overlap_len
    }

    #[inline]
    pub(crate) fn rev_start_offset(self) -> usize {
        self.r2_start_offset
    }

    #[inline]
    pub(crate) fn rev_end_offset(self) -> usize {
        self.r2_start_offset + self.overlap_len - 1
    }

    #[inline]
    pub(crate) fn reverse_range(self) -> Range<usize> {
        self.r2_start_offset..self.r2_start_offset + self.overlap_len
    }
}

#[derive(Debug, Clone)]
pub struct PairOverlap<'a> {
    slices: OrientedPairSlices<'a>,
    bounds: OverlapBounds,
}

impl<'a> PairOverlap<'a> {
    #[must_use]
    pub fn len(&self) -> usize {
        self.bounds.overlap_len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn forward_start_offset(&self) -> usize {
        self.bounds.fwd_start_offset()
    }

    #[must_use]
    pub fn forward_end_offset(&self) -> usize {
        self.bounds.fwd_end_offset()
    }

    #[must_use]
    pub fn reverse_start_offset(&self) -> usize {
        self.bounds.rev_start_offset()
    }

    #[must_use]
    pub fn reverse_end_offset(&self) -> usize {
        self.bounds.rev_end_offset()
    }

    #[must_use]
    pub fn forward_sequence(&self) -> &[u8] {
        &self.slices.forward_sequence()[self.bounds.forward_range()]
    }

    #[must_use]
    pub fn forward_qualities(&self) -> &[u8] {
        &self.slices.forward_quality_score_bytes()[self.bounds.forward_range()]
    }

    #[must_use]
    pub fn reverse_sequence(&self) -> &[u8] {
        &self.slices.reverse_sequence_rc()[self.bounds.reverse_range()]
    }

    #[must_use]
    pub fn reverse_qualities(&self) -> &[u8] {
        &self.slices.reverse_quality_score_bytes_rc()[self.bounds.reverse_range()]
    }

    #[must_use]
    pub fn overlap_windows(&self) -> (&[u8], &[u8]) {
        (self.forward_sequence(), self.reverse_sequence())
    }

    #[must_use]
    pub fn overlap_quality_windows(&self) -> (&[u8], &[u8]) {
        (self.forward_qualities(), self.reverse_qualities())
    }

    #[inline]
    pub(crate) fn id(&self) -> &str {
        self.slices.pair_id()
    }

    #[inline]
    pub(crate) fn forward_mate_sequence(&self) -> &[u8] {
        self.slices.forward_sequence()
    }

    #[inline]
    pub(crate) fn forward_mate_qualities(&self) -> &[u8] {
        self.slices.forward_quality_score_bytes()
    }

    #[inline]
    pub(crate) fn reverse_mate_sequence_rc(&self) -> &[u8] {
        self.slices.reverse_sequence_rc()
    }

    #[inline]
    pub(crate) fn reverse_mate_qualities_rc(&self) -> &[u8] {
        self.slices.reverse_quality_score_bytes_rc()
    }

    #[inline]
    pub(crate) fn oriented_slices(&self) -> &OrientedPairSlices<'a> {
        &self.slices
    }

    #[inline]
    pub(crate) fn bounds(&self) -> OverlapBounds {
        self.bounds
    }

    pub(crate) fn from_oriented_slices(
        slices: OrientedPairSlices<'a>,
        bounds: OverlapBounds,
    ) -> Result<Self> {
        slices.validate_overlap_bounds(bounds)?;
        Ok(Self { slices, bounds })
    }

    fn from_span(slices: OrientedPairSlices<'a>, span: OverlapSpan) -> Result<Self> {
        Self::from_oriented_slices(slices, span.bounds())
    }
}

#[cfg(test)]
mod tests {
    use crate::{Error, SequenceRead, errors::OverlapError};
    use proptest::{collection::vec, prelude::*};

    use super::*;

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

            let overlap_diff_max = compute_overlap_diff_max(overlap_len, params);

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

            let overlap_diff_max = compute_overlap_diff_max(overlap_len, params);

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
        let slices = mates.to_oriented_slices();
        OverlapFinder::new(params).scan_for_overlap_span_both(&slices)
    }

    fn scan_bounds_from_start(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Result<Option<OverlapSpan>> {
        let slices = mates.to_oriented_slices();
        let finder = OverlapFinder::new(params);
        let (read1, read2) = slices.sequences();
        finder.scan_from_start(read1, read2)
    }

    fn scan_bounds_from_end(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Result<Option<OverlapSpan>> {
        let slices = mates.to_oriented_slices();
        let finder = OverlapFinder::new(params);
        let (read1, read2) = slices.sequences();
        finder.scan_from_end(read1, read2)
    }

    #[test]
    fn test_scan_matches_oracle_from_start_simple_case() {
        // reverse-complement of this read is itself; easy to reason about expected overlap.
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
        // r2 reverse-complements to TTTTACGTACGT, so a reverse-direction scan should find
        // the r1 overlap against r2_rc starting after an initial 4-base offset.
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

        let params = OverlapParams::default().with_settings(
            2, // diff max
            3, // min overlap
            0.2, 3, // min comparisons
        );

        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");
        let overlap = scan_bounds_from_start(&mates, &params)
            .expect("from-start scanner should not error when checking canonical bounds");

        assert!(overlap.is_some());
        let span = overlap.expect("expected overlap in canonical from-start bounds fixture");
        assert_eq!(span.bounds().fwd_start_offset(), 0);
        assert_eq!(span.r1_end_inclusive(), 3);
        assert_eq!(span.bounds().rev_start_offset(), 0);
        assert_eq!(span.r2_end_inclusive(), 3);
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
        assert_eq!(span.r1_end_inclusive(), 7);
        assert_eq!(span.bounds().rev_start_offset(), 0);
        assert_eq!(span.r2_end_inclusive(), 7);
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
                r2: "TCGTACGT", // revcomp is ACGTACGA, one mismatch vs r1
                params: OverlapParams::default().with_settings(1, 7, 0.2, 8),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "diff_above_threshold_rejected",
                r1: "ACGTACGT",
                r2: "TCGTACGT", // same pair, but no mismatches allowed
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

        let got = mates.overlap(&params);
        assert!(matches!(
            got,
            Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
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
                prop_assert!(bounds.fwd_start_offset() <= hit.r1_end_inclusive());
                prop_assert!(bounds.rev_start_offset() <= hit.r2_end_inclusive());
                prop_assert!(hit.r1_end_inclusive() < r1.len());
                prop_assert!(hit.r2_end_inclusive() < r2.len());
                prop_assert_eq!(hit.r1_end_inclusive() - bounds.fwd_start_offset() + 1, hit.overlap_len());
                prop_assert_eq!(hit.r2_end_inclusive() - bounds.rev_start_offset() + 1, hit.overlap_len());
            }

        }
    }
}
