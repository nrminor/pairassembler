use std::ops::Range;

use crate::{Result, errors::OverlapError};

/// Canonical overlap span representation used by overlap scanning internals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OverlapSpan {
    bounds: OverlapBounds,
    diff: usize,
}

impl OverlapSpan {
    pub(crate) fn new(
        bounds: OverlapBounds,
        diff: usize,
        r1_len: usize,
        r2_len: usize,
    ) -> Result<Self> {
        let overlap_len = bounds.overlap_len();
        if overlap_len == 0 {
            return Err(OverlapError::InvalidOverlapLength {
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
            return Err(OverlapError::IndexOutOfBounds {
                read: "r1",
                index: r1_end - 1,
                length: r1_len,
            }
            .into());
        }

        if r2_end > r2_len {
            return Err(OverlapError::IndexOutOfBounds {
                read: "r2",
                index: r2_end - 1,
                length: r2_len,
            }
            .into());
        }

        if diff > overlap_len {
            return Err(OverlapError::InvalidOverlapLength {
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
    pub(crate) fn bounds(self) -> OverlapBounds {
        self.bounds
    }

    #[inline]
    pub(crate) fn overlap_len(self) -> usize {
        self.bounds.overlap_len()
    }

    #[inline]
    pub(crate) fn diff(self) -> usize {
        self.diff
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
