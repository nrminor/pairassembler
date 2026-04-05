//! The `validation` module handles finding and validating potential overlaps between mated
//! pairs of Illumina reads. Within the intended flow of data through `libpairassembly`, validation
//! should take place after overlapping with the API in the module `overlap`.
//!
//! Taking inspiration from the pre-merging validation in Brian Bushnell's BBMerge utility,
//! `validation` includes the `BaseCallValidator` struct, which uses a k-mer complexity heuristic to
//! determine how many overlap-facing bases must be present before an overlap is informative enough
//! to trust. Historical BBMerge-inspired documentation often refers to this as an entropy-based
//! check, but the implemented heuristic is better described as a complexity score rather than
//! Shannon entropy.

use std::array::IntoIter;

use rayon::prelude::*;
use tracing::warn;

use crate::{
    ReadPair, Result, SequenceRead,
    errors::ValidationError::{ExcessiveObservedMismatchRate, InsufficientOverlapLength},
    overlap::PairOverlap,
};
use wide::{CmpEq, f32x8, u8x16, u8x32};

const SIMD_LANES: usize = 32;
const ERROR_SIMD_LANES: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationPreset {
    Loose,
    Normal,
    Strict,
}

#[derive(Debug, Clone, Copy)]
pub struct ValidationPolicy {
    k: usize,
    strictness: Strictness,
    min_overlap_floor: usize,
    mismatch_multiplier: f32,
    mismatch_offset: f32,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct BaseCallValidator {
    policy: ValidationPolicy,
}

#[derive(Debug, Clone)]
pub struct ValidationMetrics {
    overlap_len: usize,
    min_informative_overlap_len: usize,
    mismatch_count: usize,
    expected_overlap_error_count: f32,
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
    pub const LOOSE_STRICTNESS_COMPLEXITY: usize = 30;
    pub const NORMAL_STRICTNESS_COMPLEXITY: usize = 39;
    pub const STRICT_STRICTNESS_COMPLEXITY: usize = 44;

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
                    "Extremely large complexity score {:?} requested, which is larger than the usual maximum of 44. This will likely lead to artifactual exclusion of many valid overlaps. Use results with caution.",
                    val
                );
                Strictness::Extreme(val)
            },
            _ if val > 55 => {
                // NOTE: This may eventually be adjusted to narrow down values that users can specify.
                // Custom errors may be useful here too.
                warn!(
                    "The requested complexity score of {val} is uncharted territory; normally values between 30 and 45 are used, with 39 usually being the sweet spot. Results with this value should be regarded with suspicion."
                );
                Strictness::Other(val)
            },
            _ if val > 0 => {
                warn!(
                    "The requested complexity score of {val} is uncharted territory; normally values between 30 and 45 are used, with 39 usually being the sweet spot. Results with this value should be regarded with suspicion."
                );
                Strictness::Other(val)
            },
            _ => {
                warn!(
                    "Invalid complexity score {:?} requested. Falling back to the strictness mode 'Normal', which defaults to a complexity score of 39.",
                    val,
                );
                Strictness::Normal(Self::NORMAL_STRICTNESS_COMPLEXITY)
            },
        }
    }
}

impl Default for ValidationPolicy {
    fn default() -> Self {
        Self::from_preset(ValidationPreset::Normal)
    }
}

impl ValidationPolicy {
    #[must_use]
    pub fn from_preset(preset: ValidationPreset) -> Self {
        match preset {
            ValidationPreset::Loose => Self {
                k: 3,
                strictness: Strictness::Loose(Strictness::LOOSE_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 5,
                mismatch_multiplier: 8.0,
                mismatch_offset: 1.0,
            },
            ValidationPreset::Normal => Self {
                k: 3,
                strictness: Strictness::Normal(Strictness::NORMAL_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 8,
                mismatch_multiplier: 6.0,
                mismatch_offset: 0.75,
            },
            ValidationPreset::Strict => Self {
                k: 3,
                strictness: Strictness::Strict(Strictness::STRICT_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 11,
                mismatch_multiplier: 4.0,
                mismatch_offset: 0.5,
            },
        }
    }

    #[must_use]
    pub fn with_k(self, k: usize) -> Self {
        Self { k, ..self }
    }

    #[must_use]
    fn with_strictness(self, strictness: Strictness) -> Self {
        Self { strictness, ..self }
    }

    #[must_use]
    pub fn k(&self) -> usize {
        self.k
    }

    #[must_use]
    pub fn min_overlap_floor(&self) -> usize {
        self.min_overlap_floor
    }

    #[must_use]
    pub fn mismatch_multiplier(&self) -> f32 {
        self.mismatch_multiplier
    }

    #[must_use]
    pub fn mismatch_offset(&self) -> f32 {
        self.mismatch_offset
    }

    #[must_use]
    fn strictness(self) -> Strictness {
        self.strictness
    }
}

impl BaseCallValidator {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn from_preset(preset: ValidationPreset) -> Self {
        Self {
            policy: ValidationPolicy::from_preset(preset),
        }
    }

    #[must_use]
    pub fn with_k(self, k: usize) -> Self {
        Self {
            policy: self.policy.with_k(k),
        }
    }

    fn with_strictness(self, strictness: Strictness) -> Self {
        Self {
            policy: self.policy.with_strictness(strictness),
        }
    }

    #[must_use]
    pub fn with_min_complexity_score(self, min_complexity_score: usize) -> Self {
        let min_complexity_score = Strictness::new_from_val(min_complexity_score);
        self.with_strictness(min_complexity_score)
    }

    #[must_use]
    pub fn with_min_entropy(self, min_entropy: usize) -> Self {
        self.with_min_complexity_score(min_entropy)
    }

    pub(crate) fn measure(&self, mates: &ReadPair, overlap: &PairOverlap) -> ValidationMetrics {
        let min_informative_overlap_len = self.compute_min_informative_overlap(mates);
        let mismatch_count = overlap.count_mismatches();
        let expected_overlap_error_count = self.sum_expected_overlap_errors(overlap);

        ValidationMetrics::new(
            overlap.len(),
            min_informative_overlap_len,
            mismatch_count,
            expected_overlap_error_count,
        )
    }

    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub(crate) fn evaluate(&self, metrics: &ValidationMetrics) -> Result<()> {
        let k = self.policy.k();
        let min_complexity_score = self.policy.strictness().get();

        let min_overlap_len = metrics
            .min_informative_overlap_len()
            .max(self.policy.min_overlap_floor());

        if metrics.overlap_len() < min_overlap_len {
            return Err(InsufficientOverlapLength {
                observed_overlap_len: metrics.overlap_len(),
                min_overlap_len,
                min_entropy: min_complexity_score,
                k,
            }
            .into());
        }

        let allowed_mismatches = metrics.expected_overlap_error_count()
            * self.policy.mismatch_multiplier()
            + self.policy.mismatch_offset();
        let observed_mismatch_count = metrics.mismatch_count() as f32;

        if observed_mismatch_count > allowed_mismatches {
            return Err(ExcessiveObservedMismatchRate {
                min_entropy: min_complexity_score,
                k,
                observed_error_rate: metrics.observed_error_rate(),
                maximum_expected_error_rate: allowed_mismatches / metrics.overlap_len() as f32,
            }
            .into());
        }

        Ok(())
    }

    /// Assess whether an overlap satisfies the configured validation policy.
    ///
    /// # Errors
    ///
    /// Returns an error if the measured overlap metrics fail the configured validation policy.
    pub fn assess(&self, mates: &ReadPair, overlap: &PairOverlap) -> Result<ValidationMetrics> {
        let metrics = self.measure(mates, overlap);
        self.evaluate(&metrics)?;
        Ok(metrics)
    }

    #[must_use]
    pub fn compute_min_informative_overlap(&self, mates: &ReadPair) -> usize {
        // pull out parameters for readability
        let k = self.policy.k();
        let min_score = self.policy.strictness().get();

        // run some sanity checks
        assert!(
            min_score <= (1 << (2 * k)),
            "min_score ({min_score}) too high for k-mer size {k}"
        );
        let read1_len = mates.fwd_sequence().len();
        let read2_len = mates.rev_sequence().len();
        assert!(
            k <= read1_len && k <= read2_len,
            "k-mer size ({k}) must not exceed read lengths (r1: {read1_len}, r2: {read2_len})"
        );

        // create mutable overlap containers for each thread to avoid data races and borrow checker gotchas
        let mut read1_head_min = 0;
        let mut read2_head_min = 0;
        let mut read1_tail_min = 0;
        let mut read2_tail_min = 0;

        // compute the minimum informative overlap using k-mer complexity from both ends of each read,
        // all
        // in parallel thanks to rayon
        rayon::scope(|s| {
            s.spawn(|_| {
                read1_head_min =
                    utils::min_overlap_by_complexity_head(mates.fwd_sequence_bytes(), k, min_score);
            });
            s.spawn(|_| {
                read2_head_min =
                    utils::min_overlap_by_complexity_head(mates.rev_sequence_bytes(), k, min_score);
            });
            s.spawn(|_| {
                read1_tail_min =
                    utils::min_overlap_by_complexity_tail(mates.fwd_sequence_bytes(), k, min_score);
            });
            s.spawn(|_| {
                read2_tail_min =
                    utils::min_overlap_by_complexity_tail(mates.rev_sequence_bytes(), k, min_score);
            });
        });

        // use whichever minimum number of overlapping bases is the highest between the pair of reads
        // and preserve the sentinel behavior from the complexity scanners when the threshold is
        // never met.
        let minimum_overlap = read1_head_min
            .max(read1_tail_min)
            .max(read2_head_min)
            .max(read2_tail_min);

        let max_possible_overlap = read1_len.min(read2_len);
        minimum_overlap.min(max_possible_overlap + 1)
    }

    #[must_use]
    pub fn compute_min_overlap(&self, mates: &ReadPair) -> usize {
        self.compute_min_informative_overlap(mates)
    }

    fn sum_expected_overlap_errors(&self, overlap: &PairOverlap) -> f32 {
        // make sure the sequence and quality lengths are correct for the forward overlap window
        let fwd_seq = overlap.forward_sequence();
        let fwd_qual = overlap.forward_qualities();
        assert_eq!(fwd_seq.len(), fwd_qual.len());

        // same for the reverse overlap window
        let rev_seq = overlap.reverse_sequence();
        let rev_qual = overlap.reverse_qualities();
        assert_eq!(rev_seq.len(), rev_qual.len());

        // initialize variables for the sums of forward and reverse errors and file with parallel
        // threads from rayon
        let mut fwd_sum_errors = 0.;
        let mut rev_sum_errors = 0.;
        rayon::scope(|s| {
            s.spawn(|_| fwd_sum_errors = sum_errors_simd(fwd_seq, fwd_qual, true));
            s.spawn(|_| rev_sum_errors = sum_errors_simd(rev_seq, rev_qual, true));
        });

        // use the lower error count as a cutoff
        fwd_sum_errors.min(rev_sum_errors)
    }
}

impl ValidationMetrics {
    pub(crate) fn new(
        overlap_len: usize,
        min_informative_overlap_len: usize,
        mismatch_count: usize,
        expected_overlap_error_count: f32,
    ) -> Self {
        Self {
            overlap_len,
            min_informative_overlap_len,
            mismatch_count,
            expected_overlap_error_count,
        }
    }

    #[must_use]
    pub fn overlap_len(&self) -> usize {
        self.overlap_len
    }

    #[must_use]
    pub fn min_informative_overlap_len(&self) -> usize {
        self.min_informative_overlap_len
    }

    #[must_use]
    pub fn min_overlap_len(&self) -> usize {
        self.min_informative_overlap_len()
    }

    #[must_use]
    pub fn mismatch_count(&self) -> usize {
        self.mismatch_count
    }

    #[must_use]
    pub fn expected_overlap_error_count(&self) -> f32 {
        self.expected_overlap_error_count
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn observed_error_rate(&self) -> f32 {
        self.mismatch_count as f32 / self.overlap_len as f32
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn expected_overlap_error_rate(&self) -> f32 {
        self.expected_overlap_error_count / self.overlap_len as f32
    }
}

impl<'overlap> PairOverlap<'overlap> {
    fn count_mismatches(&self) -> usize {
        let overlap_len = self.len();
        debug_assert_eq!(
            self.forward_end_offset() + 1 - self.forward_start_offset(),
            overlap_len
        );
        debug_assert_eq!(
            self.reverse_end_offset() + 1 - self.reverse_start_offset(),
            overlap_len
        );

        let mismatch_count =
            count_mismatches_simd(self.forward_sequence(), self.reverse_sequence());

        debug_assert!(mismatch_count < overlap_len);
        mismatch_count
    }

    fn compute_error_rate(&self) -> f32 {
        let mismatch_count = self.count_mismatches() as f32;
        let overlap_len = self.len() as f32;
        mismatch_count / overlap_len
    }

    /// Validate this overlap against the provided pair and validator policy.
    ///
    /// # Errors
    ///
    /// Returns an error if the overlap is too short or exceeds the configured mismatch/error-rate
    /// policy for the provided validator.
    pub fn validate(
        self,
        mates: &'overlap ReadPair<'overlap>,
        validator: &BaseCallValidator,
    ) -> Result<ValidatedOverlap<'overlap>> {
        let metrics = validator.assess(mates, &self)?;
        let validated = ValidatedOverlap::new_unchecked(mates, self, metrics);
        Ok(validated)
    }
}

fn count_mismatches_simd(left: &[u8], right: &[u8]) -> usize {
    let compare_len = left.len().min(right.len());
    let mut mismatches = 0usize;

    let (left_chunks, left_tail) = left[..compare_len].as_chunks::<SIMD_LANES>();
    let (right_chunks, right_tail) = right[..compare_len].as_chunks::<SIMD_LANES>();

    debug_assert_eq!(left_tail.len(), right_tail.len());

    for idx in 0..left_chunks.len() {
        let left_vec = u8x32::from(left_chunks[idx]);
        let right_vec = u8x32::from(right_chunks[idx]);
        let equal_mask = left_vec.simd_eq(right_vec).to_bitmask();
        mismatches += SIMD_LANES - equal_mask.count_ones() as usize;
    }

    for idx in 0..left_tail.len() {
        mismatches += usize::from(left_tail[idx] != right_tail[idx]);
    }

    mismatches
}

fn sum_errors_simd(seq: &[u8], qual: &[u8], count_undefined: bool) -> f32 {
    debug_assert_eq!(seq.len(), qual.len());

    let mut total = 0.0f32;
    let (seq_chunks, seq_tail) = seq.as_chunks::<ERROR_SIMD_LANES>();
    let (qual_chunks, qual_tail) = qual.as_chunks::<ERROR_SIMD_LANES>();

    debug_assert_eq!(seq_tail.len(), qual_tail.len());

    for idx in 0..seq_chunks.len() {
        let probs = error_prob_chunk(seq_chunks[idx], qual_chunks[idx], count_undefined);
        total += horizontal_sum_f32x8(probs);
    }

    for idx in 0..seq_tail.len() {
        total += scalar_error_prob(seq_tail[idx], qual_tail[idx], count_undefined);
    }

    total
}

fn error_prob_chunk(
    seq: [u8; ERROR_SIMD_LANES],
    qual: [u8; ERROR_SIMD_LANES],
    count_undefined: bool,
) -> f32x8 {
    let mut probs = [0.0f32; ERROR_SIMD_LANES];

    for idx in 0..ERROR_SIMD_LANES {
        probs[idx] = scalar_error_prob(seq[idx], qual[idx], count_undefined);
    }

    f32x8::from(probs)
}

#[inline]
fn horizontal_sum_f32x8(values: f32x8) -> f32 {
    let lanes: [f32; ERROR_SIMD_LANES] = values.into();
    lanes.into_iter().sum()
}

#[inline]
fn scalar_error_prob(base: u8, qual: u8, count_undefined: bool) -> f32 {
    if is_defined_base(base) {
        utils::phred_to_error_prob(qual)
    } else if count_undefined {
        utils::phred_to_error_prob(0)
    } else {
        0.0
    }
}

#[inline]
fn is_defined_base(base: u8) -> bool {
    let lanes = [base; 16];
    let vec = u8x16::from(lanes);
    let a = vec.simd_eq(u8x16::from([b'A'; 16])).to_bitmask() != 0;
    let c = vec.simd_eq(u8x16::from([b'C'; 16])).to_bitmask() != 0;
    let g = vec.simd_eq(u8x16::from([b'G'; 16])).to_bitmask() != 0;
    let t = vec.simd_eq(u8x16::from([b'T'; 16])).to_bitmask() != 0;
    a || c || g || t
}

#[derive(Debug)]
pub struct ValidatedOverlap<'read> {
    mates: &'read ReadPair<'read>,
    overlap: PairOverlap<'read>,
    metrics: ValidationMetrics,
}

impl<'read> ValidatedOverlap<'read> {
    pub(crate) fn new_unchecked(
        mates: &'read ReadPair<'read>,
        overlap: PairOverlap<'read>,
        metrics: ValidationMetrics,
    ) -> Self {
        Self {
            mates,
            overlap,
            metrics,
        }
    }

    #[must_use]
    pub fn read_pair(&self) -> &'read ReadPair<'read> {
        self.mates
    }

    #[must_use]
    pub fn overlap(&self) -> &PairOverlap<'read> {
        &self.overlap
    }

    #[must_use]
    pub fn validation_metrics(&self) -> &ValidationMetrics {
        &self.metrics
    }

    fn try_new(overlap: PairOverlap<'read>, mates: &'read ReadPair<'read>) -> Result<Self> {
        let validator = BaseCallValidator::default();
        let validated = overlap.validate(mates, &validator)?;
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
        let pair = self.read_pair();
        let read1 = pair.forward_read();
        let read2 = pair.reverse_read();
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

    pub(super) fn min_overlap_by_complexity_head(
        bases: &[u8],
        k: usize,
        min_score: usize,
    ) -> usize {
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

    pub(super) fn min_overlap_by_complexity_tail(bases: &[u8], k: usize, minscore: usize) -> usize {
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

#[cfg(test)]
mod tests {
    use super::{utils, *};
    use crate::{Error, errors::ValidationError};

    fn perfect_pair_fixture() -> ReadPair<'static> {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAA";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let r1 = SequenceRead::new("read1", seq, qual);
        let r2 = SequenceRead::new("read1", seq, qual);
        ReadPair::from(r1, r2).expect("test fixture reads should share the same id")
    }

    fn full_length_overlap_fixture<'a>(mates: &'a ReadPair<'a>) -> PairOverlap<'a> {
        let seq = mates.fwd_sequence_bytes();
        let qual = mates.fwd_quality_bytes();

        PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq,
            qual,
            seq.to_vec(),
            qual.to_vec(),
        )
        .expect("full-overlap fixture should satisfy overlap invariants")
    }

    #[test]
    fn test_encode_kmer_rejects_invalid_bases() {
        assert!(utils::encode_kmer(b"ACN").is_none());
        assert!(utils::encode_kmer(b"XYZ").is_none());
    }

    #[test]
    fn test_complexity_min_overlap_bounds() {
        let seq = b"ACGTACGTACGTACGT";
        let k = 3;
        let min_score = 20;

        let head = utils::min_overlap_by_complexity_head(seq, k, min_score);
        let tail = utils::min_overlap_by_complexity_tail(seq, k, min_score);

        assert!(head >= k);
        assert!(head <= seq.len() + 1);
        assert!(tail >= k);
        assert!(tail <= seq.len() + 1);
    }

    #[test]
    fn test_compute_min_informative_overlap_is_bounded() {
        let r1 = SequenceRead::new("read1", "ACGTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGTACGT", "IIIIIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let validator = BaseCallValidator::new()
            .with_k(3)
            .with_min_complexity_score(39);
        let min_overlap = validator.compute_min_informative_overlap(&mates);

        assert!(min_overlap >= 1);
        assert!(min_overlap <= mates.fwd_sequence().len().min(mates.rev_sequence().len()) + 1);
    }

    #[test]
    fn test_compute_min_informative_overlap_preserves_complexity_sentinel() {
        // Low-complexity reads plus a very strict complexity score are likely to trigger the
        // internal `len + 1` sentinel in complexity scanning.
        let r1 = SequenceRead::new("read1", "AAAAAAAAAAAA", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "AAAAAAAA", "IIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let validator = BaseCallValidator::new()
            .with_k(3)
            .with_min_complexity_score(55);
        let min_overlap = validator.compute_min_informative_overlap(&mates);

        assert_eq!(
            min_overlap,
            mates.fwd_sequence().len().min(mates.rev_sequence().len()) + 1
        );
    }

    #[test]
    fn test_validate_accepts_perfect_overlap_with_loose_settings() {
        let mates = perfect_pair_fixture();
        let overlap = full_length_overlap_fixture(&mates);

        let validator = BaseCallValidator::from_preset(ValidationPreset::Loose);
        let validated = overlap.validate(&mates, &validator);
        assert!(validated.is_ok());
    }

    #[test]
    fn test_validate_accepts_perfect_overlap_with_normal_settings() {
        let mates = perfect_pair_fixture();
        let overlap = full_length_overlap_fixture(&mates);

        let validator = BaseCallValidator::from_preset(ValidationPreset::Normal);
        let validated = overlap.validate(&mates, &validator);

        assert!(validated.is_ok());
    }

    #[test]
    fn test_validate_rejects_low_complexity_perfect_overlap_under_loose_preset() {
        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = BaseCallValidator::from_preset(ValidationPreset::Loose);
        let result = overlap.validate(&mates, &validator);

        assert!(matches!(
            result,
            Err(Error::ValidationError(
                ValidationError::InsufficientOverlapLength { .. }
            ))
        ));
    }

    #[test]
    fn test_assess_retains_validation_metrics_for_successful_overlap() {
        let mates = perfect_pair_fixture();
        let overlap = full_length_overlap_fixture(&mates);
        let validator = BaseCallValidator::from_preset(ValidationPreset::Normal);

        let metrics = validator
            .assess(&mates, &overlap)
            .expect("perfect overlap should assess successfully");

        assert_eq!(metrics.overlap_len(), overlap.len());
        assert!(metrics.min_informative_overlap_len() <= metrics.overlap_len());
        assert_eq!(metrics.mismatch_count(), 0);
        assert!(metrics.observed_error_rate().abs() < f32::EPSILON);
        assert!(metrics.expected_overlap_error_count() >= 0.0);
    }

    #[test]
    fn test_validate_rejects_short_overlap_with_insufficient_length_error() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");

        let overlap = PairOverlap::try_new(
            8,
            0,
            7,
            0,
            7,
            &seq.as_bytes()[..8],
            &qual.as_bytes()[..8],
            seq.as_bytes()[..8].to_vec(),
            qual.as_bytes()[..8].to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = BaseCallValidator::from_preset(ValidationPreset::Strict);
        let result = overlap.validate(&mates, &validator);

        assert!(matches!(
            result,
            Err(Error::ValidationError(ValidationError::InsufficientOverlapLength {
                observed_overlap_len,
                min_overlap_len,
                ..
            })) if observed_overlap_len < min_overlap_len
        ));
    }

    #[test]
    fn test_validate_rejects_excessive_mismatch_rate_in_strict_mode() {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAT";
        let mismatch_seq = "TCATTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAC";

        let r1 = SequenceRead::new("read1", seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            mismatch_seq.as_bytes().to_vec(),
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = BaseCallValidator::from_preset(ValidationPreset::Strict);
        let result = overlap.validate(&mates, &validator);
        assert!(matches!(
            result,
            Err(Error::ValidationError(ValidationError::ExcessiveObservedMismatchRate {
                observed_error_rate,
                maximum_expected_error_rate,
                ..
            })) if observed_error_rate > maximum_expected_error_rate
        ));
    }

    #[test]
    fn test_loose_preset_accepts_overlap_rejected_by_strict_preset() {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAT";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mismatch_seq = "TCGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAT";
        let mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", mismatch_seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let loose_overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            qual.as_bytes(),
            mismatch_seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");
        let strict_overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            qual.as_bytes(),
            mismatch_seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let loose = BaseCallValidator::from_preset(ValidationPreset::Loose);
        let strict = BaseCallValidator::from_preset(ValidationPreset::Strict);

        assert!(loose_overlap.validate(&mates, &loose).is_ok());
        assert!(matches!(
            strict_overlap.validate(&mates, &strict),
            Err(Error::ValidationError(
                ValidationError::ExcessiveObservedMismatchRate { .. }
            ))
        ));
    }

    #[test]
    fn test_expected_overlap_error_count_increases_with_worse_overlap_qualities() {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAT";
        let high_qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let low_qual = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        let mates = ReadPair::from(
            SequenceRead::new("read1", seq, high_qual),
            SequenceRead::new("read1", seq, high_qual),
        )
        .expect("test fixture reads should share the same id");
        let high_overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            high_qual.as_bytes(),
            seq.as_bytes().to_vec(),
            high_qual.as_bytes().to_vec(),
        )
        .expect("high-quality overlap should satisfy overlap invariants");
        let low_overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            low_qual.as_bytes(),
            seq.as_bytes().to_vec(),
            low_qual.as_bytes().to_vec(),
        )
        .expect("low-quality overlap should satisfy overlap invariants");

        let validator = BaseCallValidator::from_preset(ValidationPreset::Strict);
        let high_metrics = validator.measure(&mates, &high_overlap);
        let low_metrics = validator.measure(&mates, &low_overlap);

        assert!(
            low_metrics.expected_overlap_error_count()
                > high_metrics.expected_overlap_error_count()
        );
    }

    #[test]
    fn test_simd_mismatch_count_matches_scalar_oracle() {
        let left =
            b"ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAATCGTAGCTAGTACGATCG";
        let right =
            b"ACGTTGCAGATCTGTCCTGAATCGTACGAGTCTAGCATATGCTAGTCGATCGTACATGATCGAATCGTAGCTAGTTCGATCG";

        let observed = count_mismatches_simd(left, right);
        let expected = left
            .iter()
            .zip(right.iter())
            .filter(|(l, r)| l != r)
            .count();

        assert_eq!(observed, expected);
    }

    #[test]
    fn test_simd_expected_error_sum_matches_scalar_oracle() {
        let seq =
            b"ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAATCGTAGCTAGTACGATCG";
        let qual = "I".repeat(seq.len()).into_bytes();

        let observed = sum_errors_simd(seq, &qual, true);
        let expected = utils::sum_errors(seq, &qual, true);

        assert!((observed - expected).abs() < f32::EPSILON);
    }

    #[test]
    fn test_strict_preset_never_accepts_when_loose_rejects_for_same_overlap() {
        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let loose = BaseCallValidator::from_preset(ValidationPreset::Loose);
        let strict = BaseCallValidator::from_preset(ValidationPreset::Strict);

        assert!(overlap.validate(&mates, &loose).is_err());

        let strict_overlap = PairOverlap::try_new(
            seq.len(),
            0,
            seq.len() - 1,
            0,
            seq.len() - 1,
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");

        assert!(strict_overlap.validate(&mates, &strict).is_err());
    }
}
