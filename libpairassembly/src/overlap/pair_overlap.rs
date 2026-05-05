use super::{HasOrientedPairSlices, OrientedPairSlices, OverlapBounds, OverlapSpan};
use crate::Result;

#[derive(Debug, Clone, Copy)]
pub struct PairOverlap<'pair, 'scratch> {
    slices: OrientedPairSlices<'pair, 'scratch>,
    bounds: OverlapBounds,
}

impl<'pair, 'scratch> PairOverlap<'pair, 'scratch> {
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
    pub(crate) fn oriented_slices(&self) -> &OrientedPairSlices<'pair, 'scratch> {
        &self.slices
    }

    #[inline]
    pub(crate) fn bounds(&self) -> OverlapBounds {
        self.bounds
    }

    pub(crate) fn from_oriented_slices(
        slices: OrientedPairSlices<'pair, 'scratch>,
        bounds: OverlapBounds,
    ) -> Result<Self> {
        slices.validate_overlap_bounds(bounds)?;
        Ok(Self { slices, bounds })
    }

    pub(super) fn from_span(
        slices: OrientedPairSlices<'pair, 'scratch>,
        span: OverlapSpan,
    ) -> Result<Self> {
        Self::from_oriented_slices(slices, span.bounds())
    }
}
