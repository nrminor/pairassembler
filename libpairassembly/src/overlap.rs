use crate::{ReadMates, ValidatedOverlap};

/// Soon to be deprecated
mod methods {
    mod bbmerge;
    mod fastp;
    mod vsearch;
}

/// Parameters for the overlap analysis, mostly for skipping read pairs that no need further confirmation
/// that no overlap exists
#[derive(Debug)]
pub struct OverlapParams {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// set the minimum amount of base comparisons required to determine if
    /// two reads overlap
    min_comparisons: usize,
    search_direction: SearchDirection,
}

impl Default for OverlapParams {
    fn default() -> Self {
        OverlapParams {
            overlap_diff_max: 2,
            min_overlap: 30,
            diff_percent_max: 0.2,
            min_comparisons: 50,
            search_direction: SearchDirection::FromStart,
        }
    }
}

impl OverlapParams {
    pub fn new(
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        search_direction: SearchDirection,
    ) -> Self {
        Self {
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            search_direction,
        }
    }

    pub fn with_settings(
        mut self,
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        search_orientation: SearchDirection,
    ) -> Self {
        self.overlap_diff_max = overlap_diff_max;
        self.min_overlap = min_overlap;
        self.diff_percent_max = diff_percent_max;
        self.min_comparisons = min_comparisons;
        self.search_direction = search_orientation;
        self
    }

    pub fn with_overlap_diff_max(mut self, val: usize) -> Self {
        self.overlap_diff_max = val;
        self
    }

    pub fn with_min_overlap(mut self, val: usize) -> Self {
        self.min_overlap = val;
        self
    }

    pub fn with_diff_percent_max(mut self, val: f32) -> Self {
        self.diff_percent_max = val;
        self
    }

    pub fn with_min_comparisons(mut self, val: usize) -> Self {
        self.min_comparisons = val;
        self
    }

    pub fn with_search_direction(mut self, val: SearchDirection) -> Self {
        self.search_direction = val;
        self
    }
}

impl ReadMates<'_> {
    // pub async fn try_find_overlap(
    //     &mut self,
    //     params: &mut OverlapParams,
    // ) -> color_eyre::Result<Option<MateOverlap<'_>>> {
    //     let Some(overlap_bounds) = self.overlap_both_ends(params).await? else {
    //         return Ok(None);
    //     };

    //     todo!()
    // }

    async fn overlap_both_ends(
        &mut self,
        params: &mut OverlapParams,
    ) -> color_eyre::Result<Option<RawOverlapBounds>> {
        // TODO: use spawn_blocking for CPU bound threads here

        // search for hits from the start of read 1 first
        params.search_direction = FromStart;
        let overlap_from_left = self.scan_for_overlap_bounds(params)?;

        // if the prior search failed, try from the end of read 2
        params.search_direction = FromEnd;
        let overlap_from_right = self.scan_for_overlap_bounds(params)?;

        match (overlap_from_left, overlap_from_right) {
            // if no hits were found but no errors were encountered, return nothing
            (None, None) => Ok(None),

            // if an overlap from the right end of read 2 was found, return it
            (None, Some(right_hit)) => Ok(Some(right_hit)),

            // if an overlap from the left end of read 1 was found (the most common case), return it
            (Some(left_hit), None) => Ok(Some(left_hit)),

            // if a hit from both ends was found, compare their error rates. If there's a clear winner,
            // return that one. Otherwise, return an error; something is awry.
            (Some(left_hit), Some(right_hit)) => {
                let left_error_rate = left_hit.diff / left_hit.overlap_len;
                let right_error_rate = right_hit.diff / right_hit.overlap_len;

                // use a fancy match with guards to compactly do lots of boolean logic
                match left_error_rate {
                    _ if left_error_rate < right_error_rate => Ok(Some(left_hit)),
                    _ if right_error_rate < left_error_rate => Ok(Some(right_hit)),
                    _ => Err(eyre!("")),
                }
            },
        }
    }

    #[inline]
    fn scan_for_overlap_bounds(
        &self,
        params: &OverlapParams,
    ) -> color_eyre::Result<Option<RawOverlapBounds>> {
        // pull out the reverse complements of the reads for use later
        let read1_revcomp = self.fwd_mate.reverse_complement();
        let read2_revcomp = self.rev_mate.reverse_complement();

        // Compute the lengths of the reads and assert that the read2 length
        // is the same as the length of its complement
        let read1_len = self.fwd_mate.len();
        let read2_len = self.rev_mate.len();

        // make sure read lengths make sense
        debug_assert_eq!(
            read2_len,
            read2_revcomp.len(),
            "reverse complement length mismatch"
        );
        debug_assert_eq!(read1_len, self.fwd_mate.sequence().len());
        debug_assert_eq!(read2_len, self.rev_mate.sequence().len());

        // initialize a mutable counter for  the offset index into read 1 or 2 where the overlap begins.
        // When we search from the start of read 1, this index is initially zero. When we start from the
        // end of read 2, this index is initially the final index of read 2 and is gradually decremented
        // until an acceptable overlap is found.
        let mut overlap_start = match params.search_direction {
            FromStart => 0,
            FromEnd => read2_len - 1,
        };

        // Overlapping proceeds by either moving from the start of read 1 until it has exhausted
        // potential matches, or from the end of read 2 until the same thing occurs. At each iteration,
        // the mutable offset is updated, either moving it to the right in the case of read 1, or
        // to the left in the case of read 2. At each iteration, the algorithm will attempt to expand
        // out the largest possible overlap while only permitting mismatches up to a certain maximum.
        // If the maximum is reached and an overlap is two short, the offset is incremented, and a
        // new overlap is searched for.
        while match params.search_direction {
            FromStart => overlap_start < read1_len - params.min_overlap,
            FromEnd => overlap_start > read2_len - params.min_overlap,
        } {
            // For this iteration, the overlap length we're going for is either the length of the
            // read we're working from--read 1 if we're iterating from the start or read 2 if we're
            // iterating from the end of the pair--minus the current offset, or the length of the
            // other read, whichever is smaller.
            let overlap_len = match params.search_direction {
                FromStart => {
                    // TODO: Replace with logging
                    eprintln!("Seeking an overlap from the left end of read 1...");
                    (read1_len.saturating_sub(overlap_start)).min(read2_len)
                },
                FromEnd => {
                    // TODO: Replace with logging
                    eprintln!("Seeking an overlap from the right end of read 2...");
                    (read2_len.saturating_sub(overlap_start)).min(read1_len)
                },
            };
            debug_assert!(
                overlap_len <= read1_len && overlap_len <= read2_len,
                "Overlap length too large"
            );
            debug_assert!(
                overlap_len >= params.min_overlap,
                "Loop invariant: overlap_len < min_overlap"
            );

            // break if somehow we've ended up in a situation where we're looking for an unacceptably
            // small overlap
            if overlap_len < params.min_overlap {
                break;
            }

            // There is a certain number of overlaps that is too many. That number is either the
            // overlap difference maximum preset, or the overlap length times the preset difference
            // percentage maximum. The accomodates the fact that overlaps will sometimes be long
            // enough that a higher difference maximum is necessary.
            let overlap_diff_max = params
                .overlap_diff_max
                .min((overlap_len as f32 * params.diff_percent_max) as usize);

            // Initialize mutable counters for the number of encountered differences, the number of
            // bases compared, the position in read 1, and the position in read 2.
            //
            // That's a lot of mutable state up front, btw! This is starting to look a bit like Fortran...
            let mut diff = 0;
            let mut compared = 0;
            let mut start_in_r1: usize = 0;
            let mut start_in_r2: usize = 0;
            let mut stop_in_r1: usize = 0;
            let mut stop_in_r2: usize = 0;

            // iterate through the indices of the current overlap length. In the beginning of
            // iterations, overlaps could be as much as the entire length of read 1 or read 2.
            for overlap_pos in 0..overlap_len {
                compared = overlap_pos + 1;

                // find what positions we're at for reads 1 and 2. This is a bit tricky, as we're
                // either iterating from the start of read 1 and the start of read 2, or the the
                // end of read 1 and the start of read two.
                (start_in_r1, stop_in_r1) = params.search_direction.current_r1_bounds(
                    read1_len,
                    overlap_len,
                    overlap_start,
                    overlap_pos,
                );
                (start_in_r2, stop_in_r2) = params.search_direction.current_r2_bounds(
                    read2_len,
                    overlap_len,
                    overlap_start,
                    overlap_pos,
                );

                // Run some bounds checks
                debug_assert!(start_in_r1 < read1_len, "r1 start out of bounds");
                debug_assert!(start_in_r2 < read2_len, "r2 start out of bounds");
                debug_assert!(start_in_r1 < self.fwd_mate.sequence().len());
                debug_assert!(start_in_r2 < read2_revcomp.len());

                // If there's a mismatch, add it to the tally, breaking the loop for the overlap if
                // the count of differences is already higher than the difference max defined above
                // and too few comparisons have been performed.
                if self.fwd_mate.sequence().as_bytes()[start_in_r1]
                    != read2_revcomp.as_bytes()[start_in_r2]
                {
                    diff += 1;
                    if diff > overlap_diff_max && compared < params.min_comparisons {
                        eprintln!(
                            // TODO: Replace with logging macro
                            "Breaking at {:?} because the diff {:?} is too big for a diff max of {:?} and the number of comparisons {:?} is too small for the required {:?}.",
                            overlap_start, diff, overlap_diff_max, compared, params.min_comparisons
                        );
                        break;
                    }
                }
            }

            // On the off-chance that that for-loop completed or mostly completed, check if the
            // number of differences is lesser than or equal to the maximum differences, and
            // also that enough comparisons were performed. If both are false, return an overlap
            // result with the offset, length, and difference count for this overlap.
            if diff <= overlap_diff_max && compared >= params.min_comparisons {
                debug_assert!(diff <= overlap_diff_max);
                debug_assert!(compared >= params.min_comparisons);
                debug_assert!(overlap_len >= compared, "Compared too many bases?");
                return Ok(Some(RawOverlapBounds {
                    offset: overlap_start,
                    overlap_len,
                    diff,
                    r1_start: start_in_r1,
                    r1_end: stop_in_r1,
                    r2_start: start_in_r2,
                    r2_end: stop_in_r2,
                }));
            };

            // If we didn't early-return, increment the offset if we're moving from the start
            // of the pair, or decrement if we're moving from the end.
            match params.search_direction {
                FromStart => overlap_start += 1,
                FromEnd => overlap_start -= 1,
            };
        }

        // If the whole while loop completed without early-returning, return a NoOverlap variant
        // within an Ok, as an overlap wasn't found, but no errors occurred.
        Ok(None)
    }
}

#[derive(Debug, Default)]
pub enum SearchDirection {
    #[default]
    FromStart,
    FromEnd,
}
pub use SearchDirection::*;
use color_eyre::eyre::eyre;

impl SearchDirection {
    #[inline]
    fn current_r1_bounds(
        &self,
        read1_len: usize,
        overlap_len: usize,
        overlap_start: usize,
        offset_into_overlap: usize,
    ) -> (usize, usize) {
        let (start, stop) = match self {
            // Moving from the start of read 1 is the must intuitive case; it involves "sliding" read 1
            // gradually to the left along the left end of read 2 until an overlap is found (if at all).
            // In this case, the start and stop moves gradually to the right.
            FromStart => {
                let start = overlap_start + offset_into_overlap;
                let stop = start + (overlap_len - offset_into_overlap); // end-exclusive!

                debug_assert!(start < stop);
                (start, stop)
            },

            // When moving from the end of read 1, the stop position is always the final index of the
            // read, which is read length - 1. The start is what adjusts based on how much overlap
            // with read 2 can be achieved.
            FromEnd => {
                let stop = read1_len - 1;
                let start = stop - (overlap_len - offset_into_overlap);

                debug_assert!(start < stop);
                (start, stop)
            },
        };

        assert!(start < stop);
        (start, stop)
    }

    #[inline]
    fn current_r2_bounds(
        &self,
        read2_len: usize,
        overlap_len: usize,
        overlap_start: usize,
        offset_into_overlap: usize,
    ) -> (usize, usize) {
        let (start, stop) = match self {
            // When we move from the start of read one to find an overlap, the position in read 2
            // that's currently being compared is simply the index of the based being assessed in
            // the overlap, as starting from read 1's start also means starting from read 2's start.
            // This makes it possible to find cases where entire reads overlap, though in practice this
            // may virtually never happen.
            FromStart => {
                let start = offset_into_overlap;
                let stop = start + (overlap_len - offset_into_overlap) - 1; // end-exclusive!

                debug_assert!(start < stop);
                (start, stop)
            },

            // When moving from the end of read 1, the start position in read 2
            FromEnd => {
                // recall that when we're searching from the end of read 2, the overlap_start is,
                // confusingly, an index at or close to the end of read 2. In that case, the offset
                // will be to the left of the near-end index.
                let start = read2_len.saturating_sub(overlap_len - overlap_start);
                let stop = start + (overlap_len - offset_into_overlap) - 1;

                debug_assert!(start < stop);
                (start, stop)
            },
        };

        assert!(start < stop);
        (start, stop)
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct RawOverlapBounds {
    offset: usize,
    overlap_len: usize,
    diff: usize,
    r1_start: usize,
    r1_end: usize,
    r2_start: usize,
    r2_end: usize,
}

#[derive(Debug)]
pub struct MateOverlap<'a> {
    pub overlap_len: usize,
    pub r1_start_offset: usize,
    pub r1_end_offset: usize,
    pub r2_start_offset: usize,
    pub r2_end_offset: usize,
    pub r1_seq_view: &'a [u8],
    pub r1_qual_view: &'a [u8],
    pub r2_seq_view: &'a [u8],
    pub r2_qual_view: &'a [u8],
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use crate::Read;

    use super::*;

    #[test]
    fn test_overlap_from_start_correct_bounds() {
        let r1 = Read::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = Read::new("read1", "ACGT", "IIII");

        let mut params = OverlapParams::default().with_settings(
            2, // diff max
            3, // min overlap
            0.2,
            3, // min comparisons
            SearchDirection::FromStart,
        );

        let overlap = ReadMates {
            fwd_mate: r1,
            rev_mate: r2,
        }
        .scan_for_overlap_bounds(&params)
        .unwrap();

        assert!(overlap.is_some());
        let bounds = overlap.unwrap();
        assert_eq!(bounds.r1_start, 4);
        assert_eq!(bounds.r1_end, 8);
        assert_eq!(bounds.r2_start, 0);
        assert_eq!(bounds.r2_end, 4);
    }

    #[test]
    fn test_overlap_from_end_correct_bounds() {
        let r1 = Read::new("read1", "TTTTACGT", "IIIIIIII");
        let r2 = Read::new("read1", "ACGTAAAA", "IIIIIIII");

        let mut params =
            OverlapParams::default().with_settings(1, 4, 0.1, 4, SearchDirection::FromEnd);

        let overlap = ReadMates {
            fwd_mate: r1,
            rev_mate: r2,
        }
        .scan_for_overlap_bounds(&params)
        .unwrap();

        assert!(overlap.is_some());
        let bounds = overlap.unwrap();
        assert_eq!(bounds.r1_start, 4);
        assert_eq!(bounds.r1_end, 8);
        assert_eq!(bounds.r2_start, 0);
        assert_eq!(bounds.r2_end, 4);
    }

    #[test]
    fn test_no_overlap_detected() {
        let r1 = Read::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = Read::new("read1", "TTTT", "IIII");

        let mut params =
            OverlapParams::default().with_settings(1, 4, 0.1, 4, SearchDirection::FromStart);

        let overlap = ReadMates {
            fwd_mate: r1,
            rev_mate: r2,
        }
        .scan_for_overlap_bounds(&params)
        .unwrap();

        assert!(overlap.is_none());
    }

    #[test]
    fn test_from_start_r1_bounds() {
        let orientation = FromStart;
        let (start, stop) = orientation.current_r1_bounds(100, 30, 10, 5);
        assert_eq!(start, 15);
        assert_eq!(stop, 40);
    }

    #[test]
    fn test_from_end_r1_bounds() {
        let orientation = FromEnd;
        let (start, stop) = orientation.current_r1_bounds(100, 30, 10, 5);
        assert_eq!(stop, 99);
        assert_eq!(start, 99 + 1 - (30 - 5)); // should be 75
        assert_eq!(start, 75);
    }

    #[test]
    fn test_from_start_r2_bounds() {
        let orientation = FromStart;
        let (start, stop) = orientation.current_r2_bounds(100, 30, 10, 5);
        assert_eq!(start, 5);
        assert_eq!(stop, 30);
    }

    #[test]
    fn test_from_end_r2_bounds() {
        let orientation = FromEnd;
        let (start, stop) = orientation.current_r2_bounds(100, 30, 10, 5);
        // start = 100 - (30 + 10) = 60
        assert_eq!(start, 60);
        assert_eq!(stop, 60 + (30 - 5)); // 85
        assert_eq!(stop, 85);
    }

    #[test]
    fn test_r2_bounds_do_not_underflow() {
        let orientation = FromEnd;
        let (start, stop) = orientation.current_r2_bounds(50, 40, 20, 0);
        assert!(start <= stop); // ensures saturating_sub handled it
    }
}
