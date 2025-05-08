#![allow(clippy::pedantic, clippy::perf)]
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

//! The `validation` module handles finding and validating potential overlaps between mated
//! pairs of Illumina reads. Within the intended flow of data through `libpairassembly`, validation
//! should take place after overlapping with the API in the module `overlap`.
//!
//! Taking inspiration from the pre-merging validation in Brian Bushnell's BBMerge utility,
//! `validation` includes the `BaseCallValidator` struct, which provides parameters for using the
//! informational entropy of a given read's sequence to determine how many bases should overlap
//! in order for that overlap to be trustworthy. It also provides the ability to exclude a potential
//! overlap based on the rate of mismatches (BBMerge's so-called "ratio mode") or based on the
//! raw number of mismatches (BBMerge's "flat mode"). For most use cases, though especially for
//! long reads, filtering by error rate (the number of mismatches divided by the length of the
//! overlap) is recommended and is the default behavior in the `pairasm` CLI frontend to
//! `libpairassembly`.

use color_eyre::eyre::eyre;
use rayon::prelude::*;

use crate::{Read, ReadMates, overlap::MateOverlap};

#[derive(Debug, Clone, Copy)]
pub struct BaseCallValidator {
    k: usize,
    strictness: Strictness,
}

#[derive(Debug, Clone, Copy, Default)]
enum Strictness {
    Loose,
    #[default]
    Normal,
    Strict,
    Other(usize),
}

impl Strictness {
    const LOOSE: usize = 30;
    const NORMAL: usize = 39;
    const STRICT: usize = 44;

    fn get(&self) -> usize {
        match self {
            Strictness::Loose => Self::LOOSE,
            Strictness::Normal => Self::NORMAL,
            Strictness::Strict => Self::STRICT,

            // Warning! We got a move here!
            Strictness::Other(val) => *val,
        }
    }

    fn new_from_u8(val: usize) -> Self {
        match val {
            30 => Strictness::Loose,
            39 => Strictness::Normal,
            44 => Strictness::Strict,
            _ if val > 0 => Strictness::Other(val),
            _ => Strictness::Normal,
        }
    }
}

impl Default for BaseCallValidator {
    fn default() -> Self {
        let k = 3;
        let min_entropy = Strictness::default();
        Self {
            k,
            strictness: min_entropy,
        }
    }
}

impl BaseCallValidator {
    fn new() -> Self {
        Self::default()
    }

    fn with_k(self, k: usize) -> Self {
        Self { k, ..self }
    }

    fn with_strictness(self, strictness: Strictness) -> Self {
        Self { strictness, ..self }
    }

    fn with_min_entropy(self, min_entropy: usize) -> Self {
        let min_entropy = match min_entropy {
            30 => Strictness::Loose,
            39 => Strictness::Normal,
            44 => Strictness::Strict,
            _ if min_entropy > 0 => Strictness::Other(min_entropy),
            _ => Strictness::Normal,
        };
        Self {
            strictness: min_entropy,
            ..self
        }
    }

    fn compute_min_overlap(&self, mates: &ReadMates) -> usize {
        // pull out parameters for readability
        let k = self.k;
        let min_score = self.strictness.get();

        // run some sanity checks
        assert!(
            min_score <= (1 << (2 * k)),
            "min_score ({min_score}) too high for k-mer size {k}"
        );
        let read1_len = mates.fwd_mate.len();
        let read2_len = mates.rev_mate.len();
        assert!(
            k <= read1_len && k <= read2_len,
            "k-mer size ({k}) must not exceed read lengths (r1: {read1_len}, r2: {read2_len})"
        );

        // create mutable overlap containers for each thread
        let mut read1_head_min = 0;
        let mut read2_head_min = 0;
        let mut read1_tail_min = 0;
        let mut read2_tail_min = 0;

        // compute minimum overlap using information entropy from both ends of each read in the pair, all
        // in parallel thanks to rayon
        rayon::scope(|s| {
            s.spawn(|_| {
                read1_head_min =
                    min_overlap_by_entropy_head(mates.fwd_mate.sequence().as_bytes(), k, min_score)
            });
            s.spawn(|_| {
                read2_head_min =
                    min_overlap_by_entropy_head(mates.rev_mate.sequence().as_bytes(), k, min_score)
            });
            s.spawn(|_| {
                read1_tail_min =
                    min_overlap_by_entropy_tail(mates.fwd_mate.sequence().as_bytes(), k, min_score)
            });
            s.spawn(|_| {
                read2_tail_min =
                    min_overlap_by_entropy_tail(mates.rev_mate.sequence().as_bytes(), k, min_score)
            });
        });

        // use whichever minimum number of overlapping bases is the highest between the pair of reads
        let minimum_overlap = read1_head_min
            .max(read1_tail_min)
            .max(read2_head_min)
            .max(read2_tail_min);

        assert!(
            minimum_overlap <= read1_len.max(read2_len),
            "Computed overlap ({minimum_overlap}) exceeds maximum read length"
        );

        minimum_overlap
    }

    fn sum_expected_errors(&self, mates: &ReadMates) -> f32 {
        // make sure the sequence and quality lengths are correct for the forward mate
        let fwd_seq = mates.fwd_mate.sequence().as_bytes();
        let fwd_qual = mates.fwd_mate.quality_scores().as_bytes();
        assert_eq!(fwd_seq.len(), fwd_qual.len());

        // same for the reverse mate
        let rev_seq = mates.rev_mate.sequence().as_bytes();
        let rev_qual = mates.rev_mate.quality_scores().as_bytes();
        assert_eq!(rev_seq.len(), rev_qual.len());

        // initialize variables for the sums of forward and reverse errors and file with parallel
        // threads from rayon
        let mut fwd_sum_errors = 0.;
        let mut rev_sum_errors = 0.;
        rayon::scope(|s| {
            s.spawn(|_| fwd_sum_errors = sum_errors(fwd_seq, fwd_qual, true));
            s.spawn(|_| rev_sum_errors = sum_errors(rev_seq, rev_qual, true));
        });

        // use the lower error count as a cutoff
        fwd_sum_errors.min(rev_sum_errors)
    }
}

impl<'overlap> MateOverlap<'overlap> {
    fn count_mismatches(&self) -> usize {
        debug_assert!(((self.r1_end_offset + 1) - self.r1_start_offset) == self.overlap_len);
        debug_assert!(((self.r2_end_offset + 1) - self.r2_start_offset) == self.overlap_len);

        let mismatch_count = self
            .r1_seq_view
            .iter()
            .zip(self.r2_seq_view.iter())
            .filter(|(r1_base, r2_base)| r1_base != r2_base)
            .count();

        debug_assert!(mismatch_count < self.overlap_len);
        mismatch_count
    }

    fn compute_error_rate(&self) -> f32 {
        let mismatch_count = self.count_mismatches() as f32;
        let overlap_len = self.overlap_len as f32;
        mismatch_count / overlap_len
    }

    fn try_validate(
        self,
        mates: ReadMates<'overlap>,
        // strictness,        entropy_level
    ) -> color_eyre::Result<ValidatedOverlap<'overlap>> {
        // launch a new error calculator
        let validator = BaseCallValidator::new();

        // compute a naive minimum number of overlapping bases expected for the pair given
        // the information entropy present in the two reads sequneces
        let min_overlap = validator.compute_min_overlap(&mates);

        // early return cases where too little overlap was found based on the requested strictness
        // level.
        match validator.strictness {
            // If a loose strictness (which is to say low information entropy), the minimum
            // overlap is adjusted to allow for shorter overlaps in noisier reads.
            Strictness::Loose => {
                let expected_errors = validator.sum_expected_errors(&mates);
                // (1 + Tools.min(0.04f, errorRate) * 4f)
                let adjusted = 1. + expected_errors.min(0.04) * 4.;
                let min_overlap = min_overlap as f32 * adjusted;

                if (self.overlap_len as f32) < min_overlap {
                    return Err(eyre!(""));
                }
            },

            // If we're in strict mode, make sure the observed error rate does not exceed the expected
            // error rate as well
            Strictness::Strict => {
                if self.overlap_len < min_overlap {
                    // TODO
                    return Err(eyre!(""));
                };
                let expected_error_rate =
                    validator.sum_expected_errors(&mates) / (self.overlap_len as f32);
                if self.compute_error_rate() > expected_error_rate {
                    // TODO
                    return Err(eyre!(""));
                }
            },

            // Otherwise, just make sure there's enough overlap and then proceed
            _ => {
                if self.overlap_len < min_overlap {
                    // TODO
                    return Err(eyre!(""));
                }
            },
        }

        // If the pair was not early-returned above, it passes validation and can be returned
        // in the form of a new `ValidatedOverlap` instance
        let validated = ValidatedOverlap {
            mates,
            overlap: self,
        };
        Ok(validated)
    }
}

#[derive(Debug)]
pub struct ValidatedOverlap<'read> {
    pub mates: ReadMates<'read>,
    pub overlap: MateOverlap<'read>,
}

impl<'read> ValidatedOverlap<'read> {
    fn try_new(overlap: MateOverlap<'read>, mates: ReadMates<'read>) -> color_eyre::Result<Self> {
        let validated = overlap.try_validate(mates)?;
        Ok(validated)
    }

    /// Update bases and quality scores for the two mated reads separately. This method is intended
    /// for cases when users want to error-correct their reads without fully merging them into
    /// a consensus. This is one of the areas where the boundary between the modules in this crate
    /// gets fuzzy, but it's nice to have this functionality regardless. Some people really do just
    /// want shorter reads for some reason 🤷‍♂️
    pub fn correct_unmerged(&mut self) -> &mut Self {
        todo!()
    }

    /// Method to be called on reads to extract them back out of the validation and merging process,
    /// or to pull reads out after error-correction but before merging.
    pub fn flatten_pair(self) -> [Read<'read>; 2] {
        todo!()
    }
}

use helpers::*;
mod helpers {
    use rustc_hash::FxHashSet;

    pub(super) fn encode_kmer(kmer: &[u8]) -> Option<u64> {
        let mut code = 0u64;
        for &b in kmer {
            code <<= 2;
            code |= match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None, // Invalid base
            };
        }
        Some(code)
    }

    pub(super) fn min_overlap_by_entropy_head(bases: &[u8], k: usize, min_score: usize) -> usize {
        let mut seen_once = FxHashSet::default();
        let mut seen_twice = FxHashSet::default();

        let mut singleton_count = 0;
        let mut doubleton_count = 0;

        for (i, kmer) in bases.windows(k).enumerate() {
            // let kmer = &bases[i..i + k];
            let Some(code) = encode_kmer(kmer) else {
                continue;
            };

            if !seen_once.contains(&code) && !seen_twice.contains(&code) {
                seen_once.insert(code);
                singleton_count += 1;
            } else if seen_once.contains(&code) {
                seen_once.remove(&code);
                seen_twice.insert(code);
                doubleton_count += 1;
            }

            if singleton_count * 4 + doubleton_count >= min_score {
                return i + k;
            }
        }

        bases.len() + 1
    }

    pub(super) fn min_overlap_by_entropy_tail(bases: &[u8], k: usize, minscore: usize) -> usize {
        let mut seen_once = FxHashSet::default();
        let mut seen_twice = FxHashSet::default();

        let mut singleton_count = 0;
        let mut doubleton_count = 0;

        for (i, kmer) in bases.windows(k).enumerate().rev() {
            let Some(code) = encode_kmer(kmer) else {
                continue;
            };

            if !seen_once.contains(&code) && !seen_twice.contains(&code) {
                seen_once.insert(code);
                singleton_count += 1;
            } else if seen_once.contains(&code) {
                seen_once.remove(&code);
                seen_twice.insert(code);
                doubleton_count += 1;
            }

            if singleton_count * 4 + doubleton_count >= minscore {
                return i + k;
            }
        }

        bases.len() + 1
    }

    pub(super) fn phred_to_error_prob(phred: u8) -> f32 {
        10f32.powf(-(phred as f32) / 10.0)
    }

    pub(super) fn sum_errors(seq: &[u8], qual: &[u8], count_undefined: bool) -> f32 {
        seq.iter()
            .zip(qual.iter())
            .filter_map(|(base, qual)| match base {
                b'A' | b'C' | b'G' | b'T' => Some(phred_to_error_prob(*qual)),
                _ if count_undefined => Some(phred_to_error_prob(0)),
                _ => None,
            })
            .sum()
    }
}
