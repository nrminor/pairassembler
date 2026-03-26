use crate::{
    ReadPair, Result,
    errors::OverlapError::{IndexOutOfBounds, OverlapTie, ReverseComplementLengthMismatch},
};
use tracing::info;

/// Parameters for the overlap analysis, mostly for skipping read pairs that no need further confirmation
/// that no overlap exists
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// set the minimum amount of base comparisons required to determine if two reads overlap
    min_comparisons: usize,
    search_direction: SearchDirection,
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
            search_direction: SearchDirection::FromStart,
            tie_policy: TiePolicy::PreferFromStart,
        }
    }
}

impl OverlapParams {
    // TODO: Replace this with an actual builder pattern
    #[must_use]
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
        search_orientation: SearchDirection,
    ) -> Self {
        self.overlap_diff_max = overlap_diff_max;
        self.min_overlap = min_overlap;
        self.diff_percent_max = diff_percent_max;
        self.min_comparisons = min_comparisons;
        self.search_direction = search_orientation;
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
    fn with_search_direction(mut self, val: SearchDirection) -> Self {
        self.search_direction = val;
        self
    }

    #[must_use]
    pub fn search_from_start(mut self) -> Self {
        self.search_direction = FromStart;
        self
    }

    #[must_use]
    pub fn search_from_end(mut self) -> Self {
        self.search_direction = FromEnd;
        self
    }

    #[must_use]
    pub fn with_tie_policy(mut self, tie_policy: TiePolicy) -> Self {
        self.tie_policy = tie_policy;
        self
    }
}

impl ReadPair<'_> {
    pub fn overlap(&self, params: &OverlapParams) -> Result<Option<MateOverlap<'_>>> {
        // search for overlaps at both ends, unwrapping the raw bounds of a winning overlap if found
        let Some(overlap_bounds) = self.overlap_both_ends(params)? else {
            return Ok(None);
        };

        // extract the necessary fields of the bounds with a let pattern match
        let RawOverlapBounds {
            offset: _,
            overlap_len,
            diff: _,
            r1_start,
            r1_end,
            r2_start,
            r2_end,
        } = overlap_bounds;

        // Grab slices into the forward and reverse mates
        let fwd_seq_bytes = self.fwd_mate.sequence().as_bytes();
        let fwd_qual_bytes = self.fwd_mate.quality_scores().as_bytes();

        // Reverse complement was computed in overlap logic and stored in memory, so recalculate here.
        // TODO: refactor so REVC's are computed once and there are fewer allocatons.
        let rev_seq_rc = self.rev_mate.reverse_complement(); // produces a Vec<u8>
        let rev_qual_bytes = {
            let mut initial_vec = self.rev_mate.quality_scores().as_bytes().to_vec();
            initial_vec.reverse();
            initial_vec
        };

        // Sanity check that we can slice cleanly
        if r1_end >= fwd_qual_bytes.len() {
            return Err(IndexOutOfBounds {
                read: "fwd_mate",
                index: r1_end,
                length: fwd_seq_bytes.len(),
            }
            .into());
        }
        if r2_end >= rev_qual_bytes.len() {
            return Err(IndexOutOfBounds {
                read: "rev_mate",
                index: r2_start,
                length: rev_qual_bytes.len(),
            }
            .into());
        }

        // Create and return the rich overlap struct
        let overlap = MateOverlap {
            overlap_len,
            r1_start_offset: r1_start,
            r1_end_offset: r1_end,
            r2_start_offset: r2_start,
            r2_end_offset: r2_end,
            r1_seq_view: &fwd_seq_bytes[r1_start..=r1_end],
            r1_qual_view: &fwd_qual_bytes[r1_start..=r1_end],
            // TODO: the current logic means reverse complements and slices thereof are allocated
            // a total of three times, which should be unnecessary. A refactor should reduce this
            // to at most one copy.
            r2_seq_view: rev_seq_rc[r2_start..=r2_end].to_vec(),
            r2_qual_view: rev_qual_bytes[r2_start..=r2_end].to_vec(),
        };

        Ok(Some(overlap))
    }

    fn overlap_both_ends(&self, params: &OverlapParams) -> Result<Option<RawOverlapBounds>> {
        // Initialize params for searches from both ends
        let from_start_params = OverlapParams {
            search_direction: FromStart,
            ..*params
        };
        let from_end_params = OverlapParams {
            search_direction: FromEnd,
            ..from_start_params
        };

        // allocate memory for potential overlaps from both ends and then fill it with two parallel
        // searches. We'll use temporary mutability and a rayon scope for this.
        let (overlap_from_left, overlap_from_right) = {
            let mut overlap_from_left: Result<Option<RawOverlapBounds>> = Ok(None);
            let mut overlap_from_right: Result<Option<RawOverlapBounds>> = Ok(None);

            // run the searches in parallel with a rayon scope
            rayon::scope(|s| {
                s.spawn(|_| overlap_from_left = self.scan_for_overlap_bounds(&from_start_params));
                s.spawn(|_| overlap_from_right = self.scan_for_overlap_bounds(&from_end_params));
            });

            // If overlapping didn't encounter any errors, unwrap the optional overlap and give it
            // a henceforth immutable binding
            (overlap_from_left?, overlap_from_right?)
        };

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
                let left_error_rate = left_hit.diff as f32 / left_hit.overlap_len as f32;
                let right_error_rate = right_hit.diff as f32 / right_hit.overlap_len as f32;

                // use a fancy match with guards to compactly do lots of boolean logic
                match left_error_rate {
                    _ if left_error_rate < right_error_rate => Ok(Some(left_hit)),
                    _ if right_error_rate < left_error_rate => Ok(Some(right_hit)),
                    _ => match params.tie_policy {
                        TiePolicy::Reject => Err(OverlapTie(left_error_rate).into()),
                        TiePolicy::PreferFromStart => Ok(Some(left_hit)),
                        TiePolicy::PreferFromEnd => Ok(Some(right_hit)),
                    },
                }
            },
        }
    }

    // TODO: decide if there's a good way to expose single-end searching with this and the following
    // functions in the public API
    fn overlap_from_start(&self, params: &OverlapParams) -> Result<Option<RawOverlapBounds>> {
        let params = OverlapParams {
            search_direction: FromStart,
            ..*params
        };
        match self.scan_for_overlap_bounds(&params)? {
            // if no hits were found but no errors were encountered, return nothing
            None => Ok(None),

            // if an overlap from the left end of read 1 was found (the most common case), return it
            Some(left_hit) => Ok(Some(left_hit)),
        }
    }

    fn overlap_from_end(&self, params: &OverlapParams) -> Result<Option<RawOverlapBounds>> {
        let params = OverlapParams {
            search_direction: FromEnd,
            ..*params
        };
        match self.scan_for_overlap_bounds(&params)? {
            // if no hits were found but no errors were encountered, return nothing
            None => Ok(None),

            // if an overlap from the left end of read 1 was found (the most common case), return it
            Some(left_hit) => Ok(Some(left_hit)),
        }
    }

    // TODO: It is cuckoo how big this function is
    fn scan_for_overlap_bounds(&self, params: &OverlapParams) -> Result<Option<RawOverlapBounds>> {
        // pull out the reverse complements of the reads for use later
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
        if read2_len != read2_revcomp.len() {
            return Err(ReverseComplementLengthMismatch {
                original: read2_len,
                revcomp: read2_revcomp.len(),
            }
            .into());
        }
        debug_assert_eq!(read1_len, self.fwd_mate.sequence().len());
        debug_assert_eq!(read2_len, self.rev_mate.sequence().len());

        // initialize a mutable counter for the offset where overlap alignment begins.
        // FromStart: offset into read1 (read2 always starts at 0).
        // FromEnd: offset into read2 reverse-complement (read1 always starts at 0).
        let mut overlap_start = match params.search_direction {
            FromStart => 0,
            FromEnd => 0,
        };

        let from_start_upper = read1_len.saturating_sub(params.min_overlap);
        let from_end_upper = read2_len.saturating_sub(params.min_overlap);

        // Overlapping proceeds by either moving from the start of read 1 until it has exhausted
        // potential matches, or from the end of read 2 until the same thing occurs. At each iteration,
        // the mutable offset is updated, either moving it to the right in the case of read 1, or
        // to the left in the case of read 2. At each iteration, the algorithm will attempt to expand
        // out the largest possible overlap while only permitting mismatches up to a certain maximum.
        // If the maximum is reached and an overlap is too short, the offset is incremented, and a
        // new overlap is searched for.
        while match params.search_direction {
            FromStart => overlap_start < from_start_upper,
            FromEnd => overlap_start < from_end_upper,
        } {
            // For this iteration, the overlap length we're going for is either the length of the
            // read we're working from--read 1 if we're iterating from the start or read 2 if we're
            // iterating from the end of the pair--minus the current offset, or the length of the
            // other read, whichever is smaller.
            let overlap_len = match params.search_direction {
                FromStart => {
                    info!("Seeking an overlap from the left end of read 1...");
                    (read1_len.saturating_sub(overlap_start)).min(read2_len)
                },
                FromEnd => {
                    info!("Seeking an overlap from the right end of read 2...");
                    (read2_len.saturating_sub(overlap_start)).min(read1_len)
                },
            };
            utils::validate_overlap_len(overlap_len, read1_len, read2_len, params.min_overlap)?;

            // break if somehow we've ended up in a situation where we're looking for an unacceptably
            // small overlap
            if overlap_len < params.min_overlap {
                break;
            }

            // There is a certain number of overlaps that is too many. That number is either the
            // overlap difference maximum preset, or the overlap length times the preset difference
            // percentage maximum. The accommodates the fact that overlaps will sometimes be long
            // enough that a higher difference maximum is necessary.
            let overlap_diff_max = params
                .overlap_diff_max
                .min((overlap_len as f32 * params.diff_percent_max) as usize);

            // Phew, okay, we've made it this far. We're now ready to try and expand out the overlap,
            // as by this point in the code, we know a) which direction we're moving, b) which offset
            // into the reads we're starting at, c) how long of an overlap we're attempting to find,
            // and d) how many mismatches we're willing to tolerate.
            //
            // Next, initialize mutable counters for the number of encountered differences, the number of
            // bases compared, the position in read 1, and the position in read 2.
            //
            // That's a lot of mutable state up front, btw! This is starting to look a bit like Fortran...(shudders)
            let mut diff = 0;
            let mut compared = 0;
            let mut start_in_r1: usize = 0;
            let mut start_in_r2: usize = 0;

            // iterate through each base in the potential overlap of length `overlap_len`.
            //
            // TODO:
            //
            // Replace the for loop between the dividers with a SIMD implmentation that chunks the
            // overlap into 128-bit vectors of 16 bytes so that diffs computations can become wide.
            // This implementation should the ability to short circuit, but it'll try and do so after
            // 16 comparisons are run simulatenously.
            // ------------------------------------------------------------------------------------
            for overlap_pos in 0..overlap_len {
                compared = overlap_pos + 1;

                // find what positions we're at for reads 1 and 2. This is a bit tricky, as we're
                // either iterating from the left/5' end of read 1 and the left/5' end of read 2's
                // reverse complement ~OR~ the right end/3' end of read 1 and the left/5' end of read
                // 2's reverse complement. We use a method on the current search direction to compute
                // the bounds accordingly for each read. The logic here is tricky, so this is a good
                // place to look for bugs!
                (start_in_r1, _) = params.search_direction.current_r1_bounds(
                    read1_len,
                    overlap_len,
                    overlap_start,
                    overlap_pos,
                );
                (start_in_r2, _) = params.search_direction.current_r2_bounds(
                    read2_len,
                    overlap_len,
                    overlap_start,
                    overlap_pos,
                );

                // Run bounds checks
                utils::check_bounds(start_in_r1, start_in_r2, read1_len, read2_len)?;

                // If there's a mismatch, add it to the tally, breaking the loop for the overlap if
                // the count of differences is already higher than the difference max defined above
                // and too few comparisons have been performed. Otherwise, let the loop continue
                // to the next position in the potetial overlap...and the next position..and the next
                // position...
                if self.fwd_mate.sequence().as_bytes()[start_in_r1] != read2_revcomp[start_in_r2] {
                    diff += 1;
                    if diff > overlap_diff_max && compared < params.min_comparisons {
                        info!(
                            "Breaking at {:?} because the diff {:?} is too big for a diff max of {:?} and the number of comparisons {:?} is too small for the required {:?}.",
                            overlap_start, diff, overlap_diff_max, compared, params.min_comparisons
                        );
                        break;
                    }
                }
            }
            // ------------------------------------------------------------------------------------

            // On the off-chance that the for-loop completed or mostly completed, check if the
            // number of differences is lesser than or equal to the maximum differences, and
            // also that enough comparisons were performed. If both are true, return an overlap
            // result with the offset, length, and difference count for this overlap.
            if diff <= overlap_diff_max && compared >= params.min_comparisons {
                debug_assert!(diff <= overlap_diff_max);
                debug_assert!(compared >= params.min_comparisons);
                debug_assert!(overlap_len >= compared, "Compared too many bases?");

                let (r1_start, r1_end, r2_start, r2_end) = match params.search_direction {
                    FromStart => (
                        overlap_start,
                        overlap_start + overlap_len - 1,
                        0,
                        overlap_len - 1,
                    ),
                    FromEnd => (
                        0,
                        overlap_len - 1,
                        overlap_start,
                        overlap_start + overlap_len - 1,
                    ),
                };

                return Ok(Some(RawOverlapBounds {
                    offset: overlap_start,
                    overlap_len,
                    diff,
                    r1_start,
                    r1_end,
                    r2_start,
                    r2_end,
                }));
            }

            // If we haven't early-return by now, increment the offset if we're moving from the start
            // of the pair, or decrement if we're moving from the end.
            match params.search_direction {
                FromStart => overlap_start += 1,
                FromEnd => overlap_start += 1,
            }
        }

        // If the whole while loop completed without early-returning, return a good ol' Ok None, as
        // an overlap wasn't found, but no errors occurred. Sometime we do everything right and still
        // turn up empty-handed.
        Ok(None)
    }
}

/// `SearchDirection` is a classic "fancy bool" that uses the type system to encode more specific
/// information and enforce exhaustiveness--code to be read, not just written. Like `Result` and
/// `Option` from the standard library, its variants are themselves made public so they can be
/// initialized as literals directly. Each instance of this type is just a byte, which means
/// referencing it is more expensive in terms of CPU cycles than making copies and passing them
/// by value, so we throw on a derive for `Copy` along with the usual suspects. `Copy` needs
/// `Clone`, so that's in the derive too.
#[derive(Debug, Default, Clone, Copy)]
pub enum SearchDirection {
    #[default]
    FromStart,
    FromEnd,
}
pub use SearchDirection::*;

impl SearchDirection {
    #[inline]
    fn current_r1_bounds(
        self,
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
                let stop = overlap_start + overlap_len - 1;

                debug_assert!(start <= stop);
                (start, stop)
            },

            // When moving from the end of read 1, the stop position is always the final index of the
            // read, which is read length - 1. The start is what adjusts based on how much overlap
            // with read 2 can be achieved.
            FromEnd => {
                let stop = read1_len - 1;
                let start = (read1_len - overlap_len) + offset_into_overlap;

                debug_assert!(start <= stop);
                (start, stop)
            },
        };

        assert!(start <= stop);
        (start, stop)
    }

    #[inline]
    fn current_r2_bounds(
        self,
        _read2_len: usize,
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
                let stop = overlap_len - 1;

                debug_assert!(start <= stop);
                (start, stop)
            },

            // When moving from the "end" mode, read1 starts at index 0 and we slide read2's
            // reverse complement from left to right.
            FromEnd => {
                let start = overlap_start + offset_into_overlap;
                let stop = overlap_start + overlap_len - 1;

                debug_assert!(start <= stop);
                (start, stop)
            },
        };

        assert!(start <= stop);
        (start, stop)
    }
}

/// Short-lived intermediate representation of an overlap's bounds to be used before bundling
/// with views into the overlapped reads' sequences
#[derive(Debug, PartialEq, Eq)]
struct RawOverlapBounds {
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
    pub r2_seq_view: Vec<u8>,
    pub r2_qual_view: Vec<u8>,
}

mod utils {
    use crate::{Result, errors::OverlapError};

    /// Check whether the given offsets are valid for the current sequences and return custom errors
    pub(super) fn check_bounds(
        start_in_r1: usize,
        start_in_r2: usize,
        read1_len: usize,
        read2_len: usize,
    ) -> Result<()> {
        debug_assert!(start_in_r1 < read1_len, "r1 start out of bounds");
        debug_assert!(start_in_r2 < read2_len, "r2 start out of bounds");

        if start_in_r1 >= read1_len {
            return Err(OverlapError::IndexOutOfBounds {
                read: "r1",
                index: start_in_r1,
                length: read1_len,
            }
            .into());
        }
        if start_in_r2 >= read2_len {
            return Err(OverlapError::IndexOutOfBounds {
                read: "r2",
                index: start_in_r2,
                length: read2_len,
            }
            .into());
        }
        Ok(())
    }

    /// Check whether the computed offset length is valid. If not, it may lead to index-out-of-bound
    /// errors when viewing into the compared read mates.
    pub(super) fn validate_overlap_len(
        overlap_len: usize,
        read1_len: usize,
        read2_len: usize,
        min_required: usize,
    ) -> Result<()> {
        debug_assert!(
            overlap_len <= read1_len && overlap_len <= read2_len,
            "Overlap length too large"
        );
        debug_assert!(
            overlap_len >= min_required,
            "Loop invariant: overlap_len < min_overlap"
        );

        if overlap_len > read1_len || overlap_len > read2_len || overlap_len < min_required {
            return Err(OverlapError::InvalidOverlapLength {
                computed: overlap_len,
                read1_len,
                read2_len,
                min_required,
            }
            .into());
        }

        Ok(())
    }
}
// use utils::*;

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use crate::SequenceRead;
    use proptest::prelude::*;

    use super::*;

    /// Test-only reference implementation that mirrors fastp's simple two-loop no-gap overlap
    /// search. This is intentionally explicit and non-clever so it can serve as a behavioral
    /// oracle while refactoring `scan_for_overlap_bounds`.
    fn oracle_scan_no_gap(
        mates: &ReadPair<'_>,
        params: &OverlapParams,
    ) -> Option<RawOverlapBounds> {
        let r1 = mates.fwd_mate.sequence().as_bytes();
        let r2_rc = mates.rev_mate.reverse_complement();
        let r2 = r2_rc.as_slice();

        let len1 = r1.len();
        let len2 = r2.len();

        match params.search_direction {
            FromStart => {
                let upper = len1.saturating_sub(params.min_overlap);

                for offset in 0..upper {
                    let overlap_len = (len1 - offset).min(len2);
                    if overlap_len < params.min_overlap {
                        break;
                    }

                    let overlap_diff_max = params
                        .overlap_diff_max
                        .min((overlap_len as f32 * params.diff_percent_max) as usize);

                    let mut diff = 0;
                    let mut compared = 0;

                    for i in 0..overlap_len {
                        compared = i + 1;

                        if r1[offset + i] != r2[i] {
                            diff += 1;
                            if diff > overlap_diff_max && compared < params.min_comparisons {
                                break;
                            }
                        }
                    }

                    if diff <= overlap_diff_max && compared >= params.min_comparisons {
                        return Some(RawOverlapBounds {
                            offset,
                            overlap_len,
                            diff,
                            r1_start: offset,
                            r1_end: offset + overlap_len - 1,
                            r2_start: 0,
                            r2_end: overlap_len - 1,
                        });
                    }
                }

                None
            },
            FromEnd => {
                // fastp reverse loop: offset goes 0, -1, -2, ...
                let upper = len2.saturating_sub(params.min_overlap);

                for k in 0..upper {
                    let overlap_len = len1.min(len2 - k);
                    if overlap_len < params.min_overlap {
                        break;
                    }

                    let overlap_diff_max = params
                        .overlap_diff_max
                        .min((overlap_len as f32 * params.diff_percent_max) as usize);

                    let mut diff = 0;
                    let mut compared = 0;

                    for i in 0..overlap_len {
                        compared = i + 1;

                        if r1[i] != r2[k + i] {
                            diff += 1;
                            if diff > overlap_diff_max && compared < params.min_comparisons {
                                break;
                            }
                        }
                    }

                    if diff <= overlap_diff_max && compared >= params.min_comparisons {
                        return Some(RawOverlapBounds {
                            offset: k,
                            overlap_len,
                            diff,
                            r1_start: 0,
                            r1_end: overlap_len - 1,
                            r2_start: k,
                            r2_end: k + overlap_len - 1,
                        });
                    }
                }

                None
            },
        }
    }

    #[test]
    fn test_scan_matches_oracle_from_start_simple_case() {
        // reverse-complement of this read is itself; easy to reason about expected overlap.
        let r1 = SequenceRead::new("read1", "TTTTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let mates = ReadPair {
            fwd_mate: r1,
            rev_mate: r2,
        };

        let params =
            OverlapParams::default().with_settings(2, 4, 0.2, 4, SearchDirection::FromStart);

        let expected = oracle_scan_no_gap(&mates, &params);
        let observed = mates.scan_for_overlap_bounds(&params).unwrap();

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_scan_matches_oracle_from_end_simple_case() {
        // r2 reverse-complements to TTTTACGTACGT, so a reverse-direction scan should find
        // the r1 overlap against r2_rc starting after an initial 4-base offset.
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGTAAAA", "IIIIIIIIIIII");
        let mates = ReadPair {
            fwd_mate: r1,
            rev_mate: r2,
        };

        let params = OverlapParams::default().with_settings(2, 4, 0.2, 4, SearchDirection::FromEnd);

        let expected = oracle_scan_no_gap(&mates, &params);
        let observed = mates.scan_for_overlap_bounds(&params).unwrap();

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_overlap_from_start_correct_bounds() {
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGT", "IIII");

        let params = OverlapParams::default().with_settings(
            2, // diff max
            3, // min overlap
            0.2,
            3, // min comparisons
            SearchDirection::FromStart,
        );

        let overlap = ReadPair {
            fwd_mate: r1,
            rev_mate: r2,
        }
        .scan_for_overlap_bounds(&params)
        .unwrap();

        assert!(overlap.is_some());
        let bounds = overlap.unwrap();
        assert_eq!(bounds.r1_start, 0);
        assert_eq!(bounds.r1_end, 3);
        assert_eq!(bounds.r2_start, 0);
        assert_eq!(bounds.r2_end, 3);
    }

    #[test]
    fn test_overlap_from_end_correct_bounds() {
        let r1 = SequenceRead::new("read1", "TTTTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTAAAA", "IIIIIIII");

        let params = OverlapParams::default().with_settings(1, 4, 0.1, 4, SearchDirection::FromEnd);

        let overlap = ReadPair {
            fwd_mate: r1,
            rev_mate: r2,
        }
        .scan_for_overlap_bounds(&params)
        .unwrap();

        assert!(overlap.is_some());
        let bounds = overlap.unwrap();
        assert_eq!(bounds.r1_start, 0);
        assert_eq!(bounds.r1_end, 7);
        assert_eq!(bounds.r2_start, 0);
        assert_eq!(bounds.r2_end, 7);
    }

    #[test]
    fn test_no_overlap_detected() {
        let r1 = SequenceRead::new("read1", "ACGTACGT", "IIIIIIII");
        let r2 = SequenceRead::new("read1", "TTTT", "IIII");

        let params =
            OverlapParams::default().with_settings(1, 4, 0.1, 4, SearchDirection::FromStart);

        let overlap = ReadPair {
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
        assert_eq!(stop, 39);
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
        assert_eq!(stop, 29);
    }

    #[test]
    fn test_from_end_r2_bounds() {
        let orientation = FromEnd;
        let (start, stop) = orientation.current_r2_bounds(100, 30, 10, 5);
        assert_eq!(start, 15);
        assert_eq!(stop, 39);
    }

    #[test]
    fn test_r2_bounds_do_not_underflow() {
        let orientation = FromEnd;
        let (start, stop) = orientation.current_r2_bounds(50, 40, 10, 0);
        assert!(start <= stop);
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
                name: "exact_min_overlap_at_boundary_is_not_scanned",
                r1: "GGGGACGT",
                r2: "ACGT",
                params: OverlapParams::default().with_settings(
                    0,
                    4,
                    0.0,
                    4,
                    SearchDirection::FromStart,
                ),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "below_min_overlap_rejected",
                r1: "GGGGACGT",
                r2: "ACGTA",
                params: OverlapParams::default().with_settings(
                    0,
                    5,
                    0.0,
                    4,
                    SearchDirection::FromStart,
                ),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "diff_equal_threshold_accepted",
                r1: "ACGTACGT",
                r2: "TCGTACGT", // revcomp is ACGTACGA, one mismatch vs r1
                params: OverlapParams::default().with_settings(
                    1,
                    7,
                    0.2,
                    8,
                    SearchDirection::FromStart,
                ),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "diff_above_threshold_rejected",
                r1: "ACGTACGT",
                r2: "TCGTACGT", // same pair, but no mismatches allowed
                params: OverlapParams::default().with_settings(
                    0,
                    7,
                    0.0,
                    8,
                    SearchDirection::FromStart,
                ),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "short_perfect_but_below_min_comparisons_rejected",
                r1: "ACGT",
                r2: "ACGT",
                params: OverlapParams::default().with_settings(
                    0,
                    4,
                    0.0,
                    5,
                    SearchDirection::FromStart,
                ),
                expect_overlap: false,
            },
            OverlapFixture {
                name: "asymmetric_lengths_from_end_detected",
                r1: "ACGTACGT",
                r2: "ACGTACGTAAAA",
                params: OverlapParams::default().with_settings(
                    2,
                    4,
                    0.2,
                    4,
                    SearchDirection::FromEnd,
                ),
                expect_overlap: true,
            },
            OverlapFixture {
                name: "low_complexity_no_overlap",
                r1: "AAAAAAAAAAAA",
                r2: "CCCCCCCCCCCC",
                params: OverlapParams::default().with_settings(
                    0,
                    6,
                    0.0,
                    6,
                    SearchDirection::FromStart,
                ),
                expect_overlap: false,
            },
        ];

        for fixture in fixtures {
            let q1 = "I".repeat(fixture.r1.len());
            let q2 = "I".repeat(fixture.r2.len());
            let r1 = SequenceRead::new("read1", fixture.r1, &q1);
            let r2 = SequenceRead::new("read1", fixture.r2, &q2);
            let mates = ReadPair {
                fwd_mate: r1,
                rev_mate: r2,
            };

            let expected = oracle_scan_no_gap(&mates, &fixture.params);
            let observed = mates
                .scan_for_overlap_bounds(&fixture.params)
                .unwrap_or_else(|err| {
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

    fn dna_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        proptest::collection::vec(
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
            dir in prop_oneof![Just(SearchDirection::FromStart), Just(SearchDirection::FromEnd)],
        ) {
            let max_possible_overlap = r1.len().min(r2.len());
            let min_overlap = min_overlap_raw.min(max_possible_overlap).max(1);
            let min_comparisons = min_comparisons_raw.min(max_possible_overlap).max(1);

            let params = OverlapParams::new(
                overlap_diff_max,
                min_overlap,
                diff_percent_max,
                min_comparisons,
                dir,
            );

            let q1 = "I".repeat(r1.len());
            let q2 = "I".repeat(r2.len());
            let mates = ReadPair {
                fwd_mate: SequenceRead::new("read1", &r1, &q1),
                rev_mate: SequenceRead::new("read1", &r2, &q2),
            };

            let observed = mates.scan_for_overlap_bounds(&params);
            prop_assert!(observed.is_ok(), "scanner returned unexpected error: {observed:?}");

            let observed = observed.unwrap();
            let expected = oracle_scan_no_gap(&mates, &params);
            prop_assert_eq!(&observed, &expected);

            if let Some(hit) = observed {
                prop_assert!(hit.r1_start <= hit.r1_end);
                prop_assert!(hit.r2_start <= hit.r2_end);
                prop_assert!(hit.r1_end < r1.len());
                prop_assert!(hit.r2_end < r2.len());
                prop_assert_eq!(hit.r1_end - hit.r1_start + 1, hit.overlap_len);
                prop_assert_eq!(hit.r2_end - hit.r2_start + 1, hit.overlap_len);
            }
        }
    }
}
