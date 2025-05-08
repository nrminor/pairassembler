#![allow(clippy::pedantic, clippy::perf)]
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

/*!
This module was initially a literal refactor of the source files `overlapanalysis.cpp`
and `overlapanalysis.h` from the vsearch. It was then refactored into more idiomatic
Rust to take better advantage of its type system.
*/

use color_eyre::eyre::{Result, eyre};
use std::borrow::Borrow;

use crate::{prelude::*, utils::reverse_complement};

#[derive(Debug, Default)]
pub enum OverlapResult {
    #[default]
    NoOverlap,
    Overlap {
        offset: usize,
        overlap_len: usize,
        diff: usize,
    },
}
use OverlapResult::*;

#[derive(Debug, Default)]
pub enum SearchOrientation {
    #[default]
    FromStart,
    FromEnd,
}
pub use SearchOrientation::*;

impl SearchOrientation {
    #[inline]
    fn loop_condition(
        &self,
        offset: usize,
        read1_len: usize,
        read2_len: usize,
        min_overlap: usize,
    ) -> bool {
        match self {
            SearchOrientation::FromStart => offset < read1_len - min_overlap,
            SearchOrientation::FromEnd => offset > read2_len - min_overlap,
        }
    }

    #[inline]
    fn current_r1_idx(&self, read1_len: usize, offset: usize, j: usize) -> usize {
        match self {
            FromStart => offset + j,
            FromEnd => read1_len - 1,
        }
    }

    #[inline]
    fn current_r2_idx(&self, read2_len: usize, offset: usize, j: usize) -> usize {
        match self {
            FromStart => j,
            FromEnd => offset - 1 - j,
        }
    }
}

pub struct OverlapAnalysis<'read> {
    overlap_diff_max: usize,
    min_overlap: usize,
    diff_percent_max: f32,
    /// set the minimum amount of base comparisons required to determine if
    /// two reads overlap
    min_comparisons: usize,
    search_orientation: SearchOrientation,
    read1: Option<&'read str>,
    read2: Option<&'read str>,
}

impl Default for OverlapAnalysis<'_> {
    fn default() -> Self {
        OverlapAnalysis {
            overlap_diff_max: 2,
            min_overlap: 30,
            diff_percent_max: 0.2,
            min_comparisons: 50,
            search_orientation: SearchOrientation::FromStart,
            read1: None,
            read2: None,
        }
    }
}

impl<'read> OverlapAnalysis<'read> {
    pub fn new(
        read1: &'read str,
        read2: &'read str,
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        search_orientation: SearchOrientation,
    ) -> Self {
        Self {
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            search_orientation,
            read1: Some(read1),
            read2: Some(read2),
        }
    }

    pub fn new_settings(
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        search_orientation: SearchOrientation,
    ) -> Self {
        Self {
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            search_orientation,
            ..OverlapAnalysis::default()
        }
    }

    pub fn with_settings(
        mut self,
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        search_orientation: SearchOrientation,
    ) -> Self {
        self.overlap_diff_max = overlap_diff_max;
        self.min_overlap = min_overlap;
        self.diff_percent_max = diff_percent_max;
        self.min_comparisons = min_comparisons;
        self.search_orientation = search_orientation;
        self
    }

    pub fn settings(mut self) -> Self {
        self.read1 = None;
        self.read2 = None;
        self
    }

    pub fn with_reads(mut self, read1: &'read str, read2: &'read str) -> Self {
        self.read1 = Some(read1);
        self.read2 = Some(read2);
        self
    }

    #[inline]
    fn find_overlap(&self) -> Result<OverlapResult> {
        // Make sure the user provided an input read1
        let Some(read1) = self.read1 else {
            return Err(eyre!(
                "No read1 has been set, so no overlap will be found. Use `.with_reads()` on the overlapping struct to specify reads."
            ));
        };

        // Make sure the user provided an input read2
        let Some(read2) = self.read2 else {
            return Err(eyre!(
                "No read2 has been set, so no overlap will be found. Use `.with_reads()` on the overlapping struct to specify reads."
            ));
        };

        // compute the reverse complement of read 2
        let read2_revcomp = reverse_complement(read2);

        // Compute the lengths of the reads and assert that the read2 length
        // is the same as the length of its complement
        let read1_len = read1.len();
        let read2_len = read2.len();
        assert_eq!(read2_len, read2_revcomp.len());

        // peel off a reference to the search orientation for readability
        let search_orientation = &self.search_orientation;

        // initialize a mutable counter for  the offset index into read 1 where the overlap begins
        let mut offset = 0;

        while search_orientation.loop_condition(offset, read1_len, read2_len, self.min_overlap) {
            // For this iteration, the overlap length we're going for is either the length of the
            // read we're working from--read 1 if we're iterating from the start or read 2 if we're
            // iterating from the end of the pair--minus the current offset, or the length of the
            // other read, whichever is smaller.
            let overlap_len = match search_orientation {
                FromStart => {
                    eprintln!("MOVING FROM THE START");
                    (read1_len - offset).min(read2_len)
                },
                FromEnd => {
                    eprintln!("MOVING FROM THE END");
                    (read2_len - offset).min(read1_len)
                },
            };

            // There is a certain number of overlaps that is too many. That number is either the
            // overlap difference maximum preset, or the overlap length times the preset difference
            // percentage maximum. The accomodates the fact that overlaps will sometimes be long
            // enough that a higher difference maximum is necessary.
            let overlap_diff_max = self
                .overlap_diff_max
                .min((overlap_len as f32 * self.diff_percent_max) as usize);

            // Initialize mutable counters for the number of encountered differences and the number
            // of bases compared.
            let mut diff = 0;
            let mut compared = 0;

            // iterate through the indices of the current overlap length.
            for j in 0..overlap_len {
                compared = j + 1;

                // find what positions we're at for reads 1 and 2. This is a bit tricky, as we're
                // either iterating from the start of read 1 and the start of read 2, or the the
                // end of read 1 and the start of read two.
                let position_in_r1 = search_orientation.current_r1_idx(read1_len, offset, j);
                let position_in_r2 = search_orientation.current_r2_idx(read2_len, offset, j);

                // If there's a mismatch, add it to the tally, breaking the loop for the overlap if
                // the count of differences is already higher than the difference max defined above
                // and too few comparisons have been performed.
                if read1.as_bytes()[position_in_r1] != read2_revcomp.as_bytes()[position_in_r2] {
                    diff += 1;
                    if diff > overlap_diff_max && compared < self.min_comparisons {
                        eprintln!(
                            "Breaking at {:?} because the diff {:?} is too big for a diff max of {:?} and the number of comparisons {:?} is too small for the required {:?}.",
                            offset, diff, overlap_diff_max, compared, self.min_comparisons
                        );
                        break;
                    }
                }
            }

            // On the off-chance that that for-loop completed or mostly completed, check if the
            // number of differences is lesser than or equal to the maximum differences, and
            // also that enough comparisons were performed. If both are false, return an overlap
            // result with the offset, length, and difference count for this overlap.
            if !(diff > overlap_diff_max && compared < self.min_comparisons) {
                return Ok(Overlap {
                    offset,
                    overlap_len,
                    diff,
                });
            };

            // If we didn't early-return, increment the offset if we're moving from the start
            // of the pair, or decrement if we're moving from the end.
            match self.search_orientation {
                FromStart => offset += 1,
                FromEnd => offset -= 1,
            };
        }

        // If the whole while loop completed without early-returning, return a NoOverlap variant
        // within an Ok, as an overlap wasn't found, but no errors occurred.
        Ok(NoOverlap)
    }

    pub fn overlap_both_ends(&mut self) -> Result<OverlapResult> {
        // search for hits from the start of read 1 first
        self.search_orientation = FromStart;
        let potential_overlap1 = self.find_overlap()?;
        match potential_overlap1 {
            NoOverlap => {},
            _ => return Ok(potential_overlap1),
        };

        // if the prior search failed, try from the end of read 2
        self.search_orientation = FromEnd;
        let potential_overlap2 = self.find_overlap()?;
        match potential_overlap2 {
            NoOverlap => {},
            _ => return Ok(potential_overlap2),
        };

        // if neither of those returned a hit, return NoOverlap as the result
        Ok(NoOverlap)
    }
}

impl OverlapResult {
    fn merge<'read>(self, r1: &'read str, r2: &'read str) -> Option<&'read str> {
        let Overlap {
            offset,
            overlap_len,
            diff,
        } = self
        else {
            return None;
        };

        let len1 = overlap_len + offset.max(0);
        let len2 = if offset > 0 {
            r2.len() - overlap_len
        } else {
            0
        };

        let rr2 = reverse_complement(r2);

        let merged_seq = if offset > 0 {
            &rr2[overlap_len..len2]
        } else {
            &r1[0..len1]
        };

        let read1_quals: &str = todo!();
        let read2_quals: &str = todo!();
        let merged_qual = if offset > 0 {
            &read2_quals[overlap_len..len2]
        } else {
            &read1_quals[0..len1]
        };

        let original_name: &str = todo!();
        let strand: &str = todo!();
        let merged_name = if "=".borrow() != strand {
            format!("{strand}merged_{}_{}", len1, len2)
        } else {
            format!("{original_name}merged_{}_{}", len1, len2)
        };

        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap() {
        let r1 = "CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC";
        let r2 = "ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC";
        let qual1 = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let qual2 = "#########################################################################################";

        let overlap_result = OverlapAnalysis::default()
            .with_reads(r1, r2)
            .overlap_both_ends()
            .expect("The overlap analysis encountered an error trying to find an overlap, likely because reads have not be specified.");

        match overlap_result {
            NoOverlap => panic!(
                "Expected overlap not found in the test case, indicating that either the test is incorrect or the implementation is incorrect."
            ),
            Overlap {
                offset,
                overlap_len,
                diff,
            } => {
                assert_eq!(offset, 10);
                assert_eq!(overlap_len, 79);
                assert_eq!(diff, 1);
            },
        }
    }
}
