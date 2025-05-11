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

use std::array::IntoIter;

use rayon::prelude::*;
use tracing::warn;

use crate::{
    ReadMates, Result, SequenceRead,
    errors::ValidationError::{ExcessiveObservedMismatchRate, InsufficientOverlapLength},
    overlap::MateOverlap,
};

#[derive(Debug, Clone, Copy)]
pub struct BaseCallValidator {
    k: usize,
    strictness: Strictness,
}

#[derive(Debug, Clone, Copy)]
enum Strictness {
    Loose(usize),
    Normal(usize),
    Strict(usize),
    Extreme(usize),
    Other(usize),
}

impl Default for Strictness {
    fn default() -> Self {
        Strictness::new_from_val(39)
    }
}

impl Strictness {
    fn get(&self) -> usize {
        match self {
            Strictness::Other(val)
            | Strictness::Extreme(val)
            | Strictness::Strict(val)
            | Strictness::Normal(val)
            | Strictness::Loose(val) => *val,
        }
    }

    fn new_from_val(val: usize) -> Self {
        // pattern matching in Rust is a beautiful thing
        match val {
            21..=30 => Strictness::Loose(val),
            31..=39 => Strictness::Normal(val),
            40..=44 => Strictness::Strict(val),
            45..=55 => {
                warn!(
                    "Extremely large entropy value {:?} requested, which is larger than the usual maximum of 44. This will likely lead to artifactual exclusion of many valid overlaps. Use results with caution.",
                    val
                );
                Strictness::Extreme(val)
            },
            _ if val > 55 => {
                // NOTE: This may eventually be adjusted to narrow down values that users can specify.
                // Custom errors may be useful here too.
                warn!(
                    "The requested entropy value of {val} is uncharted territory; normally values between 30 and 45 are used, with 39 usually being the sweet spot. Results with this value should be regarded with suspicion."
                );
                Strictness::Other(val)
            },
            _ if val > 0 => {
                warn!(
                    "The requested entropy value of {val} is uncharted territory; normally values between 30 and 45 are used, with 39 usually being the sweet spot. Results with this value should be regarded with suspicion."
                );
                Strictness::Other(val)
            },
            _ => {
                warn!(
                    "Invalid entropy value {:?} requested. Falling back to the strictness mode 'Normal', which defaults to an entropy value of 39.",
                    val,
                );
                Strictness::Normal(val)
            },
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
        let min_entropy = Strictness::new_from_val(min_entropy);
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

        // create mutable overlap containers for each thread to avoid data races and borrow checker gotchas
        let mut read1_head_min = 0;
        let mut read2_head_min = 0;
        let mut read1_tail_min = 0;
        let mut read2_tail_min = 0;

        // compute minimum overlap using information entropy from both ends of each read in the pair, all
        // in parallel thanks to rayon
        rayon::scope(|s| {
            s.spawn(|_| {
                read1_head_min = utils::min_overlap_by_entropy_head(
                    mates.fwd_mate.sequence().as_bytes(),
                    k,
                    min_score,
                );
            });
            s.spawn(|_| {
                read2_head_min = utils::min_overlap_by_entropy_head(
                    mates.rev_mate.sequence().as_bytes(),
                    k,
                    min_score,
                );
            });
            s.spawn(|_| {
                read1_tail_min = utils::min_overlap_by_entropy_tail(
                    mates.fwd_mate.sequence().as_bytes(),
                    k,
                    min_score,
                );
            });
            s.spawn(|_| {
                read2_tail_min = utils::min_overlap_by_entropy_tail(
                    mates.rev_mate.sequence().as_bytes(),
                    k,
                    min_score,
                );
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
            s.spawn(|_| fwd_sum_errors = utils::sum_errors(fwd_seq, fwd_qual, true));
            s.spawn(|_| rev_sum_errors = utils::sum_errors(rev_seq, rev_qual, true));
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

    pub fn try_validate(
        self,
        mates: &'overlap ReadMates<'overlap>,
        validator: &BaseCallValidator,
    ) -> Result<ValidatedOverlap<'overlap>> {
        // compute a naive minimum number of overlapping bases expected for the pair given
        // the information entropy present in the two reads sequneces
        let min_overlap_len = validator.compute_min_overlap(mates);

        // early return cases where too little overlap was found based on the requested strictness
        // level.
        match *validator {
            // If a loose strictness (which is to say low information entropy), the minimum
            // overlap is adjusted to allow for shorter overlaps in noisier reads.
            BaseCallValidator {
                k,
                strictness: Strictness::Loose(min_entropy),
            } => {
                let expected_errors = validator.sum_expected_errors(mates);
                // (1 + Tools.min(0.04f, errorRate) * 4f)
                let adjusted = 1. + expected_errors.min(0.04) * 4.;
                let min_overlap_len = min_overlap_len as f32 * adjusted;

                if (self.overlap_len as f32) < min_overlap_len {
                    return Err(InsufficientOverlapLength {
                        observed_overlap_len: self.overlap_len,
                        min_overlap_len: min_overlap_len as usize,
                        min_entropy,
                        k,
                    }
                    .into());
                }
            },

            // If we're in strict or extreme mode, make sure the observed error rate does not exceed
            // the expected error rate as well
            BaseCallValidator {
                k,
                strictness: Strictness::Strict(min_entropy) | Strictness::Extreme(min_entropy),
            } => {
                if self.overlap_len < min_overlap_len {
                    return Err(InsufficientOverlapLength {
                        observed_overlap_len: self.overlap_len,
                        min_overlap_len,
                        min_entropy,
                        k,
                    }
                    .into());
                }
                let maximum_expected_error_rate =
                    validator.sum_expected_errors(mates) / (self.overlap_len as f32);
                let observed_error_rate = self.compute_error_rate();
                if observed_error_rate > maximum_expected_error_rate {
                    return Err(ExcessiveObservedMismatchRate {
                        min_entropy,
                        k,
                        observed_error_rate,
                        maximum_expected_error_rate,
                    }
                    .into());
                }
            },

            // Otherwise, just make sure there's enough overlap and then proceed. We explicitly match
            // on the remaining variants here for futureproofing, as this means we'll get exhaustiveness
            // checking if more variants are added in the future.
            BaseCallValidator {
                k,
                strictness: Strictness::Normal(min_entropy) | Strictness::Other(min_entropy),
            } => {
                if self.overlap_len < min_overlap_len {
                    return Err(InsufficientOverlapLength {
                        observed_overlap_len: self.overlap_len,
                        min_overlap_len,
                        min_entropy,
                        k,
                    }
                    .into());
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
    pub mates: &'read ReadMates<'read>,
    pub overlap: MateOverlap<'read>,
}

impl<'read> ValidatedOverlap<'read> {
    fn try_new(overlap: MateOverlap<'read>, mates: &'read ReadMates<'read>) -> Result<Self> {
        let validator = BaseCallValidator::default();
        let validated = overlap.try_validate(mates, &validator)?;
        Ok(validated)
    }

    /// Update bases and quality scores for the two mated reads separately. This method is intended
    /// for cases when users want to error-correct their reads without fully merging them into
    /// a consensus. This is one of the areas where the boundary between the modules in this crate
    /// gets fuzzy, but it's nice to have this functionality regardless. Some people really do just
    /// want shorter reads for some reason 🤷‍♂️
    pub fn correct_unmerged(&mut self) -> &mut Self {
        unimplemented!()
    }

    /// Method to be called on reads to extract them back out of the validation and merging process,
    /// or to pull reads out after error-correction but before merging.
    #[must_use]
    pub fn extract_pair(self) -> [&'read SequenceRead<'read>; 2] {
        let read1 = &self.mates.fwd_mate;
        let read2 = &self.mates.rev_mate;
        [read1, read2]
    }
}

impl<'read> IntoIterator for ValidatedOverlap<'read> {
    type Item = &'read SequenceRead<'read>;
    type IntoIter = IntoIter<Self::Item, 2>;

    fn into_iter(self) -> Self::IntoIter {
        self.extract_pair().into_iter()
    }
}

mod utils {
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
        10f32.powf(-f32::from(phred) / 10.0)
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
