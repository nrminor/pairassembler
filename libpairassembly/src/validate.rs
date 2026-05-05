//! Overlap validation for paired sequencing reads.
//!
//! [`OverlapValidator`] uses a k-mer complexity heuristic to decide whether a discovered overlap is
//! informative enough to trust before merging or correction. Historical BBMerge-inspired
//! documentation often calls this kind of check entropy-based, but this implementation is better
//! described as a complexity score rather than Shannon entropy.

use tracing::warn;

use crate::{
    Result,
    assembler::HasPairOverlap,
    errors::ValidationError::{ExcessiveObservedMismatchRate, InsufficientOverlapLength},
    overlap::{HasOrientedPairSlices, PairOverlap},
    read::ReadPair,
};
use wide::{CmpEq, f32x8, u8x16, u8x32};

const SIMD_LANES: usize = 32;
const ERROR_SIMD_LANES: usize = 8;

/// Preset validation strictness for overlap informativeness checks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationPreset {
    /// More permissive; useful for exploratory work or shorter inserts.
    Loose,
    /// Balanced default.
    Normal,
    /// More conservative; useful when false merges are especially costly.
    Strict,
}

/// Tunable policy used by [`OverlapValidator`].
#[derive(Debug, Clone, Copy)]
pub struct ValidationPolicy {
    k: usize,
    strictness: Strictness,
    min_overlap_floor: usize,
    mismatch_multiplier: f32,
    mismatch_offset: f32,
}

/// Validates whether discovered overlap evidence is informative enough to trust.
///
/// The validator combines a k-mer complexity-derived minimum overlap length with observed and
/// expected mismatch/error rates. It is intentionally separate from overlap search: finding a
/// candidate overlap and deciding whether to trust it are different stages.
#[derive(Debug, Clone, Copy, Default)]
pub struct OverlapValidator {
    policy: ValidationPolicy,
}

/// Measurements retained from an overlap validation decision.
#[derive(Debug, Clone)]
pub struct ValidationMetrics {
    overlap_len: usize,
    min_informative_overlap_len: usize,
    mismatch_count: usize,
    expected_overlap_error_count: f64,
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
    pub const LOOSE_STRICTNESS_COMPLEXITY: usize = 24;
    pub const NORMAL_STRICTNESS_COMPLEXITY: usize = 30;
    pub const STRICT_STRICTNESS_COMPLEXITY: usize = 39;

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
    /// Build a policy from a named preset.
    ///
    /// ```rust
    /// use libpairassembly::{ValidationPolicy, ValidationPreset};
    ///
    /// let policy = ValidationPolicy::from_preset(ValidationPreset::Normal);
    /// assert_eq!(policy.k(), 3);
    /// ```
    #[must_use]
    pub fn from_preset(preset: ValidationPreset) -> Self {
        match preset {
            ValidationPreset::Loose => Self {
                k: 3,
                strictness: Strictness::Loose(Strictness::LOOSE_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 5,
                mismatch_multiplier: 10.0,
                mismatch_offset: 1.5,
            },
            ValidationPreset::Normal => Self {
                k: 3,
                strictness: Strictness::Normal(Strictness::NORMAL_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 5,
                mismatch_multiplier: 8.0,
                mismatch_offset: 1.0,
            },
            ValidationPreset::Strict => Self {
                k: 3,
                strictness: Strictness::Strict(Strictness::STRICT_STRICTNESS_COMPLEXITY),
                min_overlap_floor: 8,
                mismatch_multiplier: 6.0,
                mismatch_offset: 0.75,
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

impl OverlapValidator {
    /// Build the default overlap validator.
    ///
    /// ```rust
    /// use libpairassembly::OverlapValidator;
    ///
    /// let validator = OverlapValidator::new();
    /// let _same_default = OverlapValidator::default();
    /// # let _ = validator;
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Build an overlap validator from a named preset.
    #[must_use]
    pub fn from_preset(preset: ValidationPreset) -> Self {
        Self {
            policy: ValidationPolicy::from_preset(preset),
        }
    }

    /// Set the k-mer size used by the complexity heuristic.
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

    /// Set the minimum k-mer complexity score required by the validator.
    #[must_use]
    pub fn with_min_complexity_score(self, min_complexity_score: usize) -> Self {
        let min_complexity_score = Strictness::new_from_val(min_complexity_score);
        self.with_strictness(min_complexity_score)
    }

    /// Validate a concrete pair overlap and retain its measured validation metrics.
    ///
    /// # Errors
    ///
    /// Returns an error if the overlap is too short or exceeds the configured mismatch/error-rate
    /// policy.
    pub fn validate_overlap<'overlap, 'scratch>(
        &self,
        overlap: PairOverlap<'overlap, 'scratch>,
    ) -> Result<ValidatedOverlap<'overlap, 'scratch>> {
        let metrics = self.assess(&overlap)?;
        Ok(ValidatedOverlap::new_unchecked(overlap, metrics))
    }

    pub(crate) fn measure<T>(&self, target: &T) -> Result<ValidationMetrics>
    where
        T: HasPairOverlap + ?Sized,
    {
        target.validate_overlap_bounds()?;
        let slices = target.pair_slices()?;
        let bounds = target.overlap_bounds()?;
        let (forward_sequence, reverse_sequence) = slices.sequences();
        let (fwd_seq, rev_seq) = target.overlap_windows()?;
        let (fwd_qual, rev_qual) = target.overlap_quality_windows()?;

        let min_informative_overlap_len =
            self.compute_min_informative_overlap_for_sequences(forward_sequence, reverse_sequence);
        let mismatch_count = count_mismatches_simd(fwd_seq, rev_seq);
        let expected_overlap_error_count =
            sum_expected_overlap_errors(fwd_seq, fwd_qual, rev_seq, rev_qual);

        Ok(ValidationMetrics::new(
            bounds.overlap_len(),
            min_informative_overlap_len,
            mismatch_count,
            expected_overlap_error_count,
        ))
    }

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
                min_complexity_score,
                k,
            }
            .into());
        }

        let allowed_mismatches = metrics.expected_overlap_error_count()
            * f64::from(self.policy.mismatch_multiplier())
            + f64::from(self.policy.mismatch_offset());
        let observed_mismatch_count = usize_to_f64(metrics.mismatch_count());

        if observed_mismatch_count > allowed_mismatches {
            return Err(ExcessiveObservedMismatchRate {
                min_complexity_score,
                k,
                observed_error_rate: metrics.observed_error_rate(),
                maximum_expected_error_rate: allowed_mismatches
                    / usize_to_f64(metrics.overlap_len()),
            }
            .into());
        }

        Ok(())
    }

    /// Assess whether paired overlap slices satisfy the configured validation policy.
    ///
    /// # Errors
    ///
    /// Returns an error if the measured overlap metrics fail the configured validation policy.
    pub(crate) fn assess<T>(&self, target: &T) -> Result<ValidationMetrics>
    where
        T: HasPairOverlap + ?Sized,
    {
        let metrics = self.measure(target)?;
        self.evaluate(&metrics)?;
        Ok(metrics)
    }

    #[must_use]
    pub fn compute_min_informative_overlap(&self, mates: &ReadPair) -> usize {
        self.compute_min_informative_overlap_for_sequences(
            mates.fwd_sequence_bytes(),
            mates.rev_sequence_bytes(),
        )
    }

    #[must_use]
    fn compute_min_informative_overlap_for_sequences(&self, read1: &[u8], read2: &[u8]) -> usize {
        // pull out parameters for readability
        let k = self.policy.k();
        let min_score = self.policy.strictness().get();

        // run some sanity checks
        assert!(
            min_score <= (1 << (2 * k)),
            "min_score ({min_score}) too high for k-mer size {k}"
        );
        let read1_len = read1.len();
        let read2_len = read2.len();
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
                read1_head_min = utils::min_overlap_by_complexity_head(read1, k, min_score);
            });
            s.spawn(|_| {
                read2_head_min = utils::min_overlap_by_complexity_head(read2, k, min_score);
            });
            s.spawn(|_| {
                read1_tail_min = utils::min_overlap_by_complexity_tail(read1, k, min_score);
            });
            s.spawn(|_| {
                read2_tail_min = utils::min_overlap_by_complexity_tail(read2, k, min_score);
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
}

fn sum_expected_overlap_errors(
    fwd_seq: &[u8],
    fwd_qual: &[u8],
    rev_seq: &[u8],
    rev_qual: &[u8],
) -> f64 {
    assert_eq!(fwd_seq.len(), fwd_qual.len());
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
    f64::from(fwd_sum_errors.min(rev_sum_errors))
}

impl ValidationMetrics {
    pub(crate) fn new(
        overlap_len: usize,
        min_informative_overlap_len: usize,
        mismatch_count: usize,
        expected_overlap_error_count: f64,
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
    pub fn expected_overlap_error_count(&self) -> f64 {
        self.expected_overlap_error_count
    }

    #[must_use]
    pub fn observed_error_rate(&self) -> f64 {
        usize_to_f64(self.mismatch_count) / usize_to_f64(self.overlap_len)
    }

    #[must_use]
    pub fn expected_overlap_error_rate(&self) -> f64 {
        self.expected_overlap_error_count / usize_to_f64(self.overlap_len)
    }
}

fn usize_to_f64(value: usize) -> f64 {
    const U32_RADIX_USIZE: usize = 4_294_967_296;
    const U32_RADIX_F64: f64 = 4_294_967_296.0;

    let high = value / U32_RADIX_USIZE;
    let low = value % U32_RADIX_USIZE;

    let Ok(high) = u32::try_from(high) else {
        return f64::INFINITY;
    };
    let Ok(low) = u32::try_from(low) else {
        unreachable!("value modulo 2^32 always fits in u32")
    };

    f64::from(high) * U32_RADIX_F64 + f64::from(low)
}

impl<'overlap, 'scratch> PairOverlap<'overlap, 'scratch> {
    /// Validate this overlap against the provided validator policy.
    ///
    /// # Errors
    ///
    /// Returns an error if the overlap is too short or exceeds the configured mismatch/error-rate
    /// policy for the provided validator.
    pub fn validate(
        self,
        validator: &OverlapValidator,
    ) -> Result<ValidatedOverlap<'overlap, 'scratch>> {
        validator.validate_overlap(self)
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
pub struct ValidatedOverlap<'read, 'scratch> {
    overlap: PairOverlap<'read, 'scratch>,
    metrics: ValidationMetrics,
}

impl<'read, 'scratch> ValidatedOverlap<'read, 'scratch> {
    pub(crate) fn new_unchecked(
        overlap: PairOverlap<'read, 'scratch>,
        metrics: ValidationMetrics,
    ) -> Self {
        Self { overlap, metrics }
    }

    #[must_use]
    pub fn overlap(&self) -> &PairOverlap<'read, 'scratch> {
        &self.overlap
    }

    #[must_use]
    pub fn validation_metrics(&self) -> &ValidationMetrics {
        &self.metrics
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
                return bases.len() - i;
            }
        }

        bases.len() + 1
    }

    pub(super) fn phred_to_error_prob(phred: u8) -> f32 {
        10f32.powf(-f32::from(phred) / 10.0)
    }

    #[cfg(test)]
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
    use crate::{
        Error,
        errors::ValidationError,
        overlap::{OrientedPairSlices, OverlapBounds},
        read::SequenceRead,
    };

    fn test_overlap<'a>(
        bounds: OverlapBounds,
        fwd_seq: &'a [u8],
        fwd_qual: &'a [u8],
        rev_seq_rc: &'a [u8],
        rev_qual_rev: &'a [u8],
    ) -> Result<PairOverlap<'a, 'a>> {
        let slices = OrientedPairSlices {
            id: "read1",
            fwd_seq,
            fwd_quality_score_bytes: fwd_qual,
            rev_seq_rc,
            rev_quality_score_bytes_rc: rev_qual_rev,
        };

        PairOverlap::from_oriented_slices(slices, bounds)
    }

    fn perfect_pair_fixture() -> ReadPair<'static> {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAA";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let r1 = SequenceRead::new("read1", seq, qual);
        let r2 = SequenceRead::new("read1", seq, qual);
        ReadPair::from(r1, r2).expect("test fixture reads should share the same id")
    }

    fn full_length_overlap_fixture<'a>(mates: &'a ReadPair<'a>) -> PairOverlap<'a, 'a> {
        let seq = mates.fwd_sequence_bytes();
        let qual = mates.fwd_quality_bytes();

        test_overlap(OverlapBounds::new(seq.len(), 0, 0), seq, qual, seq, qual)
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
    fn test_tail_complexity_reports_bases_needed_from_tail() {
        let seq = b"ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAA";
        let reversed: Vec<_> = seq.iter().copied().rev().collect();
        let k = 3;

        for min_score in [30, 39, 44] {
            let from_tail = utils::min_overlap_by_complexity_tail(seq, k, min_score);
            let from_reversed_head = utils::min_overlap_by_complexity_head(&reversed, k, min_score);

            assert_eq!(
                from_tail, from_reversed_head,
                "tail complexity should report the number of bases needed from the tail, not a left-origin coordinate"
            );
        }
    }

    #[test]
    fn test_tail_complexity_requirement_is_monotonic() {
        let seq = b"ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAA";
        let k = 3;

        let loose = utils::min_overlap_by_complexity_tail(seq, k, 30);
        let normal = utils::min_overlap_by_complexity_tail(seq, k, 39);
        let strict = utils::min_overlap_by_complexity_tail(seq, k, 44);

        assert!(
            loose <= normal && normal <= strict,
            "raising the required complexity score should not reduce the required tail overlap: loose={loose}, normal={normal}, strict={strict}"
        );
    }

    #[test]
    fn test_compute_min_informative_overlap_is_bounded() {
        let r1 = SequenceRead::new("read1", "ACGTACGTACGT", "IIIIIIIIIIII");
        let r2 = SequenceRead::new("read1", "ACGTACGTACGT", "IIIIIIIIIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let validator = OverlapValidator::new()
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

        let validator = OverlapValidator::new()
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

        let validator = OverlapValidator::from_preset(ValidationPreset::Loose);
        let validated = overlap.validate(&validator);
        assert!(validated.is_ok());
    }

    #[test]
    fn test_validate_accepts_perfect_overlap_with_normal_settings() {
        let mates = perfect_pair_fixture();
        let overlap = full_length_overlap_fixture(&mates);

        let validator = OverlapValidator::from_preset(ValidationPreset::Normal);
        let validated = overlap.validate(&validator);

        assert!(validated.is_ok());
    }

    #[test]
    fn test_validate_rejects_low_complexity_perfect_overlap_under_loose_preset() {
        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let _mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes(),
            qual.as_bytes(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = OverlapValidator::from_preset(ValidationPreset::Loose);
        let result = overlap.validate(&validator);

        assert!(matches!(
            result,
            Err(Error::Validation(
                ValidationError::InsufficientOverlapLength { .. }
            ))
        ));
    }

    #[test]
    fn test_assess_retains_validation_metrics_for_successful_overlap() {
        let mates = perfect_pair_fixture();
        let overlap = full_length_overlap_fixture(&mates);
        let validator = OverlapValidator::from_preset(ValidationPreset::Normal);

        let metrics = validator
            .assess(&overlap)
            .expect("perfect overlap should assess successfully");

        assert_eq!(metrics.overlap_len(), overlap.len());
        assert!(metrics.min_informative_overlap_len() <= metrics.overlap_len());
        assert_eq!(metrics.mismatch_count(), 0);
        assert!(metrics.observed_error_rate().abs() < f64::EPSILON);
        assert!(metrics.expected_overlap_error_count() >= 0.0);
    }

    #[test]
    fn test_validate_rejects_short_overlap_with_insufficient_length_error() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let _mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");

        let overlap = test_overlap(
            OverlapBounds::new(8, 0, 0),
            &seq.as_bytes()[..8],
            &qual.as_bytes()[..8],
            &seq.as_bytes()[..8],
            &qual.as_bytes()[..8],
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = OverlapValidator::from_preset(ValidationPreset::Strict);
        let result = overlap.validate(&validator);

        assert!(matches!(
            result,
            Err(Error::Validation(ValidationError::InsufficientOverlapLength {
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

        let overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            mismatch_seq.as_bytes(),
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        )
        .expect("test overlap should satisfy overlap invariants");

        let validator = OverlapValidator::from_preset(ValidationPreset::Strict);
        let result = overlap.validate(&validator);
        assert!(matches!(
            result,
            Err(Error::Validation(ValidationError::ExcessiveObservedMismatchRate {
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
        let _mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", mismatch_seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let loose_overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            qual.as_bytes(),
            mismatch_seq.as_bytes(),
            qual.as_bytes(),
        )
        .expect("test overlap should satisfy overlap invariants");
        let strict_overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            qual.as_bytes(),
            mismatch_seq.as_bytes(),
            qual.as_bytes(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let loose = OverlapValidator::from_preset(ValidationPreset::Loose);
        let strict = OverlapValidator::from_preset(ValidationPreset::Strict);

        assert!(loose_overlap.validate(&loose).is_ok());
        assert!(matches!(
            strict_overlap.validate(&strict),
            Err(Error::Validation(
                ValidationError::ExcessiveObservedMismatchRate { .. }
            ))
        ));
    }

    #[test]
    fn test_expected_overlap_error_count_increases_with_worse_overlap_qualities() {
        let seq = "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTAT";
        let high_qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let low_qual = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        let _mates = ReadPair::from(
            SequenceRead::new("read1", seq, high_qual),
            SequenceRead::new("read1", seq, high_qual),
        )
        .expect("test fixture reads should share the same id");
        let high_overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            high_qual.as_bytes(),
            seq.as_bytes(),
            high_qual.as_bytes(),
        )
        .expect("high-quality overlap should satisfy overlap invariants");
        let low_overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            low_qual.as_bytes(),
            seq.as_bytes(),
            low_qual.as_bytes(),
        )
        .expect("low-quality overlap should satisfy overlap invariants");

        let validator = OverlapValidator::from_preset(ValidationPreset::Strict);
        let high_metrics = validator
            .measure(&high_overlap)
            .expect("high-quality overlap should measure successfully");
        let low_metrics = validator
            .measure(&low_overlap)
            .expect("low-quality overlap should measure successfully");

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
        let _mates = ReadPair::from(
            SequenceRead::new("read1", seq, qual),
            SequenceRead::new("read1", seq, qual),
        )
        .expect("test fixture reads should share the same id");
        let overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes(),
            qual.as_bytes(),
        )
        .expect("test overlap should satisfy overlap invariants");

        let loose = OverlapValidator::from_preset(ValidationPreset::Loose);
        let strict = OverlapValidator::from_preset(ValidationPreset::Strict);

        assert!(overlap.validate(&loose).is_err());

        let strict_overlap = test_overlap(
            OverlapBounds::new(seq.len(), 0, 0),
            seq.as_bytes(),
            qual.as_bytes(),
            seq.as_bytes(),
            qual.as_bytes(),
        )
        .expect("test overlap should satisfy overlap invariants");

        assert!(strict_overlap.validate(&strict).is_err());
    }
}
