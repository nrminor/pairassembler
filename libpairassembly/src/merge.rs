use crate::{
    PairOverlap, ReadPair, Result,
    assembler::{HasMergeableOverlap, IntoOwnedRecordParts},
    errors::MergeError::{
        EmptyOverlapWindow, MergeSequenceQualityLengthMismatch, MergedLengthMismatch,
        OverlapWindowLengthMismatch, ProvenanceLengthMismatch,
    },
};

/// Owned merged-layout tuple used when handing merged output to correction internals.
///
/// Layout:
/// `(id, consensus_seq, consensus_qual, left_overhang_len, fwd_overlap_seq, fwd_overlap_qual,
/// rev_overlap_seq_rc, rev_overlap_qual_rc)`.
///
/// The explicit `left_overhang_len` keeps correction aligned to the overlap window while preserving
/// overhang qualities outside the corrected region.
pub(crate) type MergeCorrectionParts = (
    String,
    Vec<u8>,
    Vec<u8>,
    usize,
    Vec<u8>,
    Vec<u8>,
    Vec<u8>,
    Vec<u8>,
);

/// Deterministic consensus read produced by the merge stage.
///
/// This is the canonical merge artifact for downstream processing.
#[derive(Debug, Clone)]
pub struct MergedRead {
    id: String,
    consensus_seq: Vec<u8>,
    consensus_qual: Vec<u8>,
    left_overhang_len: usize,
    provenance: MergeProvenance,
}

/// Overlap-only evidence retained from merging.
///
/// This payload is intentionally narrow: it records only overlap windows needed
/// for downstream correction/diagnostics and does not represent full-read
/// provenance.
#[derive(Debug, Clone)]
pub struct MergeProvenance {
    overlap_len: usize,
    fwd_overlap_seq: Vec<u8>,
    fwd_overlap_qual: Vec<u8>,
    rev_overlap_seq: Vec<u8>,
    rev_overlap_qual: Vec<u8>,
}

/// Borrowed normalized merge input projected from a pair-plus-overlap carrier.
///
/// `MergeView` is an internal normalization boundary rather than a domain object:
/// it packages the exact slices and coordinate mappings needed by the cheap merge
/// kernel without owning sequence data.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MergeView<'a> {
    id: &'a str,
    left_overhang_seq: &'a [u8],
    left_overhang_qual: &'a [u8],
    fwd_overlap_seq: &'a [u8],
    fwd_overlap_qual: &'a [u8],
    rev_raw_seq: &'a [u8],
    rev_raw_qual: &'a [u8],
    overlap_len: usize,
    rev_overlap_start_rc: usize,
    right_overhang_start_rc: usize,
}

impl<'a> MergeView<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_pair_bounds(
        pair: &'a ReadPair<'a>,
        overlap_len: usize,
        fwd_start_offset: usize,
        fwd_end_offset: usize,
        rev_start_offset: usize,
        rev_end_offset: usize,
    ) -> Result<Self> {
        if overlap_len == 0 {
            return Err(EmptyOverlapWindow.into());
        }

        let fwd_end_exclusive =
            fwd_end_offset
                .checked_add(1)
                .ok_or(OverlapWindowLengthMismatch {
                    fwd_len: 0,
                    rev_len: 0,
                })?;
        let right_overhang_start_rc =
            rev_end_offset
                .checked_add(1)
                .ok_or(OverlapWindowLengthMismatch {
                    fwd_len: 0,
                    rev_len: 0,
                })?;

        let fwd_seq = pair.fwd_sequence_bytes();
        let fwd_qual = pair.fwd_quality_bytes();
        let rev_raw_seq = pair.rev_sequence_bytes();
        let rev_raw_qual = pair.rev_quality_bytes();

        ensure_seq_qual_lengths("left", fwd_seq, fwd_qual)?;
        ensure_seq_qual_lengths("overlap_rev", rev_raw_seq, rev_raw_qual)?;

        if fwd_start_offset > fwd_end_offset
            || fwd_end_exclusive > fwd_seq.len()
            || fwd_end_exclusive > fwd_qual.len()
            || rev_start_offset > rev_end_offset
            || right_overhang_start_rc > rev_raw_seq.len()
            || right_overhang_start_rc > rev_raw_qual.len()
        {
            return Err(OverlapWindowLengthMismatch {
                fwd_len: fwd_end_exclusive.saturating_sub(fwd_start_offset),
                rev_len: right_overhang_start_rc.saturating_sub(rev_start_offset),
            }
            .into());
        }

        if (fwd_end_exclusive - fwd_start_offset) != overlap_len
            || (right_overhang_start_rc - rev_start_offset) != overlap_len
        {
            return Err(OverlapWindowLengthMismatch {
                fwd_len: fwd_end_exclusive - fwd_start_offset,
                rev_len: right_overhang_start_rc - rev_start_offset,
            }
            .into());
        }

        Ok(Self {
            id: pair.fwd_id(),
            left_overhang_seq: &fwd_seq[..fwd_start_offset],
            left_overhang_qual: &fwd_qual[..fwd_start_offset],
            fwd_overlap_seq: &fwd_seq[fwd_start_offset..fwd_end_exclusive],
            fwd_overlap_qual: &fwd_qual[fwd_start_offset..fwd_end_exclusive],
            rev_raw_seq,
            rev_raw_qual,
            overlap_len,
            rev_overlap_start_rc: rev_start_offset,
            right_overhang_start_rc,
        })
    }

    pub(crate) fn from_pair_overlap(
        pair: &'a ReadPair<'a>,
        overlap: &PairOverlap<'_>,
    ) -> Result<Self> {
        Self::from_pair_bounds(
            pair,
            overlap.len(),
            overlap.forward_start_offset(),
            overlap.forward_end_offset(),
            overlap.reverse_start_offset(),
            overlap.reverse_end_offset(),
        )
    }

    #[inline]
    pub fn id(&self) -> &str {
        self.id
    }

    #[inline]
    pub fn left_overhang_len(&self) -> usize {
        self.left_overhang_seq.len()
    }

    #[inline]
    pub fn overlap_len(&self) -> usize {
        self.overlap_len
    }

    #[inline]
    pub fn right_overhang_len(&self) -> usize {
        self.rev_raw_seq
            .len()
            .saturating_sub(self.right_overhang_start_rc)
    }

    #[inline]
    pub fn merged_len(&self) -> usize {
        self.left_overhang_len() + self.overlap_len() + self.right_overhang_len()
    }

    pub fn copy_left_overhang_seq_into(&self, out: &mut [u8]) -> usize {
        let n = self.left_overhang_seq.len().min(out.len());
        out[..n].copy_from_slice(&self.left_overhang_seq[..n]);
        n
    }

    pub fn copy_left_overhang_qual_into(&self, out: &mut [u8]) -> usize {
        let n = self.left_overhang_qual.len().min(out.len());
        out[..n].copy_from_slice(&self.left_overhang_qual[..n]);
        n
    }

    pub fn copy_fwd_overlap_seq_into(&self, out: &mut [u8]) -> usize {
        let n = self.fwd_overlap_seq.len().min(out.len());
        out[..n].copy_from_slice(&self.fwd_overlap_seq[..n]);
        n
    }

    pub fn copy_fwd_overlap_qual_into(&self, out: &mut [u8]) -> usize {
        let n = self.fwd_overlap_qual.len().min(out.len());
        out[..n].copy_from_slice(&self.fwd_overlap_qual[..n]);
        n
    }

    pub fn copy_rev_overlap_seq_rc_into(&self, out: &mut [u8]) -> usize {
        let n = self.overlap_len.min(out.len());
        for (i, slot) in out[..n].iter_mut().enumerate() {
            *slot = self.rev_overlap_base_rc_at(i);
        }
        n
    }

    pub fn copy_rev_overlap_qual_rc_into(&self, out: &mut [u8]) -> usize {
        let n = self.overlap_len.min(out.len());
        for (i, slot) in out[..n].iter_mut().enumerate() {
            *slot = self.rev_overlap_qual_rc_at(i);
        }
        n
    }

    pub fn copy_right_overhang_seq_rc_into(&self, out: &mut [u8]) -> usize {
        let n = self.right_overhang_len().min(out.len());
        for (i, slot) in out[..n].iter_mut().enumerate() {
            *slot = self.right_overhang_base_rc_at(i);
        }
        n
    }

    pub fn copy_right_overhang_qual_rc_into(&self, out: &mut [u8]) -> usize {
        let n = self.right_overhang_len().min(out.len());
        for (i, slot) in out[..n].iter_mut().enumerate() {
            *slot = self.right_overhang_qual_rc_at(i);
        }
        n
    }

    #[inline]
    pub(crate) fn left_overhang_seq(&self) -> &[u8] {
        self.left_overhang_seq
    }

    #[inline]
    pub(crate) fn left_overhang_qual(&self) -> &[u8] {
        self.left_overhang_qual
    }

    #[inline]
    pub(crate) fn fwd_overlap_seq(&self) -> &[u8] {
        self.fwd_overlap_seq
    }

    #[inline]
    pub(crate) fn fwd_overlap_qual(&self) -> &[u8] {
        self.fwd_overlap_qual
    }

    #[inline]
    pub(crate) fn overlap_pair_at(&self, i: usize) -> ((u8, u8), (u8, u8)) {
        debug_assert!(i < self.overlap_len);
        let fwd = (self.fwd_overlap_seq[i], self.fwd_overlap_qual[i]);
        let rev = (
            self.rev_overlap_base_rc_at(i),
            self.rev_overlap_qual_rc_at(i),
        );
        (fwd, rev)
    }

    #[inline]
    pub(crate) fn right_overhang_pair_at(&self, i: usize) -> (u8, u8) {
        (
            self.right_overhang_base_rc_at(i),
            self.right_overhang_qual_rc_at(i),
        )
    }

    #[inline]
    fn rev_overlap_base_rc_at(&self, i: usize) -> u8 {
        let rc_idx = self.rev_overlap_start_rc + i;
        let raw_idx = self.rc_to_raw_index(rc_idx);
        complement_base(self.rev_raw_seq[raw_idx])
    }

    #[inline]
    fn rev_overlap_qual_rc_at(&self, i: usize) -> u8 {
        let rc_idx = self.rev_overlap_start_rc + i;
        let raw_idx = self.rc_to_raw_index(rc_idx);
        self.rev_raw_qual[raw_idx]
    }

    #[inline]
    fn right_overhang_base_rc_at(&self, i: usize) -> u8 {
        let rc_idx = self.right_overhang_start_rc + i;
        let raw_idx = self.rc_to_raw_index(rc_idx);
        complement_base(self.rev_raw_seq[raw_idx])
    }

    #[inline]
    fn right_overhang_qual_rc_at(&self, i: usize) -> u8 {
        let rc_idx = self.right_overhang_start_rc + i;
        let raw_idx = self.rc_to_raw_index(rc_idx);
        self.rev_raw_qual[raw_idx]
    }

    #[inline]
    fn rc_to_raw_index(&self, rc_idx: usize) -> usize {
        debug_assert!(rc_idx < self.rev_raw_seq.len());
        self.rev_raw_seq.len() - 1 - rc_idx
    }
}

#[inline]
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

impl MergedRead {
    /// Construct a merged read without re-checking invariants.
    ///
    /// This is reserved for callers that have already validated consensus and provenance layout.
    pub(crate) fn new_unchecked(
        id: String,
        consensus_seq: Vec<u8>,
        consensus_qual: Vec<u8>,
        left_overhang_len: usize,
        provenance: MergeProvenance,
    ) -> Self {
        Self {
            id,
            consensus_seq,
            consensus_qual,
            left_overhang_len,
            provenance,
        }
    }

    /// Construct a merged read from checked consensus parts plus merge provenance.
    ///
    /// # Errors
    ///
    /// Returns an error if consensus sequence and quality lengths differ, or if the retained
    /// overlap window cannot fit within the consensus layout described by `left_overhang_len`.
    pub(crate) fn try_new(
        id: String,
        consensus_seq: Vec<u8>,
        consensus_qual: Vec<u8>,
        left_overhang_len: usize,
        provenance: MergeProvenance,
    ) -> Result<Self> {
        ensure_seq_qual_lengths("consensus", &consensus_seq, &consensus_qual)?;
        if left_overhang_len + provenance.overlap_len() > consensus_seq.len() {
            return Err(MergedLengthMismatch {
                expected: left_overhang_len + provenance.overlap_len(),
                actual: consensus_seq.len(),
            }
            .into());
        }

        Ok(Self::new_unchecked(
            id,
            consensus_seq,
            consensus_qual,
            left_overhang_len,
            provenance,
        ))
    }

    /// Borrow the merged read identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Borrow merged consensus sequence bytes.
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        self.consensus_seq.as_slice()
    }

    /// Borrow merged consensus quality bytes.
    #[must_use]
    pub fn qualities(&self) -> &[u8] {
        self.consensus_qual.as_slice()
    }

    /// Return the merged sequence length.
    #[must_use]
    pub fn len(&self) -> usize {
        self.consensus_seq.len()
    }

    /// Return whether the merged sequence is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.consensus_seq.is_empty()
    }

    /// Consume and return owned merged sequence bytes.
    #[must_use]
    pub fn sequence_owned(self) -> Vec<u8> {
        self.consensus_seq
    }

    /// Consume and return owned merged quality bytes.
    #[must_use]
    pub fn qualities_owned(self) -> Vec<u8> {
        self.consensus_qual
    }

    /// Borrow overlap evidence retained from merge.
    #[must_use]
    pub fn provenance(&self) -> &MergeProvenance {
        &self.provenance
    }

    /// Consume into the owned layout expected by correction internals.
    #[must_use]
    pub(crate) fn into_correction_parts(self) -> MergeCorrectionParts {
        let MergedRead {
            id,
            consensus_seq,
            consensus_qual,
            left_overhang_len,
            provenance,
        } = self;

        (
            id,
            consensus_seq,
            consensus_qual,
            left_overhang_len,
            provenance.fwd_overlap_seq,
            provenance.fwd_overlap_qual,
            provenance.rev_overlap_seq,
            provenance.rev_overlap_qual,
        )
    }
}

impl MergeProvenance {
    /// Construct checked overlap evidence retained from merge.
    ///
    /// # Errors
    ///
    /// Returns an error if sequence/quality lengths differ within either overlap window, or if the
    /// retained forward/reverse overlap sequence lengths disagree with `overlap_len`.
    pub(crate) fn try_new(
        overlap_len: usize,
        fwd_overlap_seq: Vec<u8>,
        fwd_overlap_qual: Vec<u8>,
        rev_overlap_seq: Vec<u8>,
        rev_overlap_qual: Vec<u8>,
    ) -> Result<Self> {
        ensure_seq_qual_lengths("overlap_fwd", &fwd_overlap_seq, &fwd_overlap_qual)?;
        ensure_seq_qual_lengths("overlap_rev", &rev_overlap_seq, &rev_overlap_qual)?;

        if fwd_overlap_seq.len() != overlap_len || rev_overlap_seq.len() != overlap_len {
            return Err(ProvenanceLengthMismatch {
                overlap_len,
                fwd_len: fwd_overlap_seq.len(),
                rev_len: rev_overlap_seq.len(),
            }
            .into());
        }

        Ok(Self {
            overlap_len,
            fwd_overlap_seq,
            fwd_overlap_qual,
            rev_overlap_seq,
            rev_overlap_qual,
        })
    }

    /// Return overlap length used by merge.
    #[must_use]
    pub fn overlap_len(&self) -> usize {
        self.overlap_len
    }

    /// Borrow forward overlap sequence bytes.
    #[must_use]
    pub fn fwd_overlap_seq(&self) -> &[u8] {
        self.fwd_overlap_seq.as_slice()
    }

    /// Borrow forward overlap quality bytes.
    #[must_use]
    pub fn fwd_overlap_qual(&self) -> &[u8] {
        self.fwd_overlap_qual.as_slice()
    }

    /// Borrow reverse overlap sequence bytes in reverse-complement orientation.
    #[must_use]
    pub fn rev_overlap_seq(&self) -> &[u8] {
        self.rev_overlap_seq.as_slice()
    }

    /// Borrow reverse overlap quality bytes in reverse-complement orientation.
    #[must_use]
    pub fn rev_overlap_qual(&self) -> &[u8] {
        self.rev_overlap_qual.as_slice()
    }
}

impl IntoOwnedRecordParts for MergedRead {
    fn into_owned_record_parts(self) -> (String, Vec<u8>, Vec<u8>) {
        (self.id, self.consensus_seq, self.consensus_qual)
    }
}

/// Merge deterministic consensus output from any carrier implementing merge capabilities.
///
/// This is the thin adapter from carrier/state types into the normalized [`MergeView`] boundary.
/// The merge kernel itself operates on `MergeView` so validation can stay upstream and the kernel
/// can remain computationally cheap.
///
/// # Errors
///
/// Returns an error if overlap windows cannot be projected into a consistent
/// merge view or if final merged-length invariants are violated.
pub(crate) fn merge_from<T>(input: &T) -> Result<MergedRead>
where
    T: HasMergeableOverlap,
{
    let view = input.merge_view()?;
    merge_kernel(view)
}

fn merge_kernel(view: MergeView<'_>) -> Result<MergedRead> {
    let overlap_len = view.overlap_len();
    debug_assert_eq!(overlap_len, view.fwd_overlap_seq().len());

    let expected_len = view.merged_len();
    let mut full_seq = Vec::with_capacity(expected_len);
    let mut full_qual = Vec::with_capacity(expected_len);

    full_seq.extend_from_slice(view.left_overhang_seq());
    full_qual.extend_from_slice(view.left_overhang_qual());

    for i in 0..overlap_len {
        let ((fb, fq), (rb, rq)) = view.overlap_pair_at(i);
        if fq >= rq {
            full_seq.push(fb);
            full_qual.push(fq);
        } else {
            full_seq.push(rb);
            full_qual.push(rq);
        }
    }

    for i in 0..view.right_overhang_len() {
        let (base, qual) = view.right_overhang_pair_at(i);
        full_seq.push(base);
        full_qual.push(qual);
    }

    if full_seq.len() != expected_len {
        return Err(MergedLengthMismatch {
            expected: expected_len,
            actual: full_seq.len(),
        }
        .into());
    }

    let mut rev_overlap_seq = vec![0u8; overlap_len];
    let mut rev_overlap_qual = vec![0u8; overlap_len];
    view.copy_rev_overlap_seq_rc_into(&mut rev_overlap_seq);
    view.copy_rev_overlap_qual_rc_into(&mut rev_overlap_qual);

    let provenance = MergeProvenance::try_new(
        overlap_len,
        view.fwd_overlap_seq().to_vec(),
        view.fwd_overlap_qual().to_vec(),
        rev_overlap_seq,
        rev_overlap_qual,
    )?;

    MergedRead::try_new(
        view.id().to_owned(),
        full_seq,
        full_qual,
        view.left_overhang_len(),
        provenance,
    )
}

#[inline]
fn ensure_seq_qual_lengths(section: &'static str, seq: &[u8], qual: &[u8]) -> Result<()> {
    if seq.len() != qual.len() {
        return Err(MergeSequenceQualityLengthMismatch {
            section,
            seq_len: seq.len(),
            qual_len: qual.len(),
        }
        .into());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{MergeView, merge_from};
    use crate::{
        Error, PairOverlap, ReadPair, SequenceRead,
        assembler::HasMergeableOverlap,
        errors::MergeError,
        validate::{ValidatedOverlap, ValidationMetrics},
    };
    use proptest::{collection::vec, prelude::*};

    fn dna_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        vec(
            prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
            min_len..=max_len,
        )
        .prop_map(|chars| chars.into_iter().collect())
    }

    fn qual_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        vec(33u8..=73u8, min_len..=max_len)
            .prop_map(|bytes| bytes.into_iter().map(char::from).collect())
    }

    #[derive(Debug, Clone)]
    struct MergeFixture {
        left_seq: String,
        overlap_fwd_seq: String,
        overlap_rev_seq: String,
        right_seq: String,
        left_qual: String,
        overlap_fwd_qual: String,
        overlap_rev_qual: String,
        right_qual: String,
    }

    prop_compose! {
        fn merge_fixture_strategy()
            (left_len in 0usize..=16, overlap_len in 4usize..=24, right_len in 0usize..=16)
            (
                left_seq in dna_string_strategy(left_len, left_len),
                overlap_fwd_seq in dna_string_strategy(overlap_len, overlap_len),
                overlap_rev_seq in dna_string_strategy(overlap_len, overlap_len),
                right_seq in dna_string_strategy(right_len, right_len),
                left_qual in qual_string_strategy(left_len, left_len),
                overlap_fwd_qual in qual_string_strategy(overlap_len, overlap_len),
                overlap_rev_qual in qual_string_strategy(overlap_len, overlap_len),
                right_qual in qual_string_strategy(right_len, right_len),
            ) -> MergeFixture
        {
            MergeFixture {
                left_seq,
                overlap_fwd_seq,
                overlap_rev_seq,
                right_seq,
                left_qual,
                overlap_fwd_qual,
                overlap_rev_qual,
                right_qual,
            }
        }
    }

    fn reverse_complement_dna(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|base| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                invalid => panic!("invalid DNA base in merge test fixture: {invalid}"),
            })
            .collect()
    }

    fn build_validated_overlap_from_fixture(fixture: &MergeFixture) -> ValidatedOverlap<'static> {
        let left_seq = fixture.left_seq.as_str();
        let overlap_fwd_seq = fixture.overlap_fwd_seq.as_str();
        let overlap_rev_seq = fixture.overlap_rev_seq.as_str();
        let right_seq = fixture.right_seq.as_str();
        let left_qual = fixture.left_qual.as_str();
        let overlap_fwd_qual = fixture.overlap_fwd_qual.as_str();
        let overlap_rev_qual = fixture.overlap_rev_qual.as_str();
        let right_qual = fixture.right_qual.as_str();

        let fwd_seq = format!("{left_seq}{overlap_fwd_seq}");
        let fwd_qual = format!("{left_qual}{overlap_fwd_qual}");

        let rev_rc_seq = format!("{overlap_rev_seq}{right_seq}");
        let rev_rc_qual = format!("{overlap_rev_qual}{right_qual}");

        let rev_seq = reverse_complement_dna(&rev_rc_seq);
        let rev_qual = rev_rc_qual.chars().rev().collect::<String>();

        let fwd_static: &'static str = Box::leak(fwd_seq.into_boxed_str());
        let fwd_qual_static: &'static str = Box::leak(fwd_qual.into_boxed_str());
        let rev_static: &'static str = Box::leak(rev_seq.into_boxed_str());
        let rev_qual_static: &'static str = Box::leak(rev_qual.into_boxed_str());

        let mates = ReadPair::from(
            SequenceRead::new("read1", fwd_static, fwd_qual_static),
            SequenceRead::new("read1", rev_static, rev_qual_static),
        )
        .expect("test fixtures should produce valid paired reads");
        let mates_ref: &'static ReadPair<'static> = Box::leak(Box::new(mates));

        let overlap_len = overlap_fwd_seq.len();
        let fwd_start = left_seq.len();
        let fwd_end = fwd_start + overlap_len - 1;
        let overlap = PairOverlap::try_new(
            overlap_len,
            fwd_start,
            fwd_end,
            0,
            overlap_len - 1,
            &mates_ref.fwd_sequence_bytes()[fwd_start..=fwd_end],
            &mates_ref.fwd_quality_bytes()[fwd_start..=fwd_end],
            overlap_rev_seq.as_bytes().to_vec(),
            overlap_rev_qual.as_bytes().to_vec(),
        )
        .expect("merge fixture overlap should satisfy overlap invariants");

        let metrics = ValidationMetrics::new(overlap_len, overlap_len, 0, 0.0);

        ValidatedOverlap::new_unchecked(mates_ref, overlap, metrics)
    }

    fn oracle_merge(fixture: &MergeFixture) -> (Vec<u8>, Vec<u8>) {
        let left_seq = fixture.left_seq.as_str();
        let overlap_fwd_seq = fixture.overlap_fwd_seq.as_str();
        let overlap_rev_seq = fixture.overlap_rev_seq.as_str();
        let right_seq = fixture.right_seq.as_str();
        let left_qual = fixture.left_qual.as_str();
        let overlap_fwd_qual = fixture.overlap_fwd_qual.as_str();
        let overlap_rev_qual = fixture.overlap_rev_qual.as_str();
        let right_qual = fixture.right_qual.as_str();

        let mut overlap_seq = Vec::with_capacity(overlap_fwd_seq.len());
        let mut overlap_qual = Vec::with_capacity(overlap_fwd_qual.len());

        for (((fwd_base, fwd_q), rev_base), rev_q) in overlap_fwd_seq
            .bytes()
            .zip(overlap_fwd_qual.bytes())
            .zip(overlap_rev_seq.bytes())
            .zip(overlap_rev_qual.bytes())
        {
            if fwd_q >= rev_q {
                overlap_seq.push(fwd_base);
                overlap_qual.push(fwd_q);
            } else {
                overlap_seq.push(rev_base);
                overlap_qual.push(rev_q);
            }
        }

        let mut full_seq = Vec::new();
        full_seq.extend_from_slice(left_seq.as_bytes());
        full_seq.extend_from_slice(&overlap_seq);
        full_seq.extend_from_slice(right_seq.as_bytes());

        let mut full_qual = Vec::new();
        full_qual.extend_from_slice(left_qual.as_bytes());
        full_qual.extend_from_slice(&overlap_qual);
        full_qual.extend_from_slice(right_qual.as_bytes());

        (full_seq, full_qual)
    }

    #[test]
    fn test_merge_perfect_full_overlap_roundtrip() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = PairOverlap::try_new(
            5,
            4,
            8,
            0,
            4,
            &mates.fwd_sequence_bytes()[4..=8],
            &mates.fwd_quality_bytes()[4..=8],
            crate::prelude::utils::reverse_complement(mates.rev_sequence()).into_bytes(),
            mates.rev_quality_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");
        let metrics = ValidationMetrics::new(5, 5, 0, 0.0);
        let validated = ValidatedOverlap::new_unchecked(&mates, overlap, metrics);

        let merged = merge_from(&validated)
            .expect("generic merge_from should merge validated overlap without bounds errors");

        assert_eq!(merged.id(), "read1");
        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.qualities(), b"IIIIIIIII");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_with_left_overhang_preserves_prefix() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = PairOverlap::try_new(
            5,
            4,
            8,
            0,
            4,
            &mates.fwd_sequence_bytes()[4..=8],
            &mates.fwd_quality_bytes()[4..=8],
            crate::prelude::utils::reverse_complement(mates.rev_sequence()).into_bytes(),
            mates.rev_quality_bytes().to_vec(),
        )
        .expect("test overlap should satisfy overlap invariants");
        let metrics = ValidationMetrics::new(5, 5, 0, 0.0);
        let validated = ValidatedOverlap::new_unchecked(&mates, overlap, metrics);

        let merged = merge_from(&validated)
            .expect("generic merge_from should merge validated overlap without bounds errors");

        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_with_right_overhang_preserves_suffix() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: "TT".into(),
            overlap_fwd_seq: "ACGT".into(),
            overlap_rev_seq: "ACGT".into(),
            right_seq: "GG".into(),
            left_qual: "II".into(),
            overlap_fwd_qual: "IIII".into(),
            overlap_rev_qual: "IIII".into(),
            right_qual: "II".into(),
        });

        let merged = merge_from(&validated)
            .expect("generic merge_from should preserve right overhang semantics");

        assert_eq!(merged.sequence(), b"TTACGTGG");
        assert_eq!(merged.qualities(), b"IIIIIIII");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_tie_on_quality_prefers_forward_base() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: String::new(),
            overlap_fwd_seq: "AAAA".into(),
            overlap_rev_seq: "TTTT".into(),
            right_seq: String::new(),
            left_qual: String::new(),
            overlap_fwd_qual: "IIII".into(),
            overlap_rev_qual: "IIII".into(),
            right_qual: String::new(),
        });

        let merged = merge_from(&validated)
            .expect("generic merge_from should remain deterministic on equal-quality overlap");

        assert_eq!(merged.sequence(), b"AAAA");
        assert_eq!(merged.qualities(), b"IIII");
    }

    #[test]
    fn test_merge_from_populates_provenance_overlap_evidence() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: "TT".into(),
            overlap_fwd_seq: "ACGT".into(),
            overlap_rev_seq: "TCGT".into(),
            right_seq: "GG".into(),
            left_qual: "II".into(),
            overlap_fwd_qual: "IIII".into(),
            overlap_rev_qual: "IIII".into(),
            right_qual: "II".into(),
        });

        let migrated = merge_from(&validated)
            .expect("generic merge_from should succeed on validated overlap fixture");

        assert_eq!(migrated.provenance().overlap_len(), 4);
        assert_eq!(migrated.provenance().fwd_overlap_seq(), b"ACGT");
        assert_eq!(migrated.provenance().fwd_overlap_qual(), b"IIII");
        assert_eq!(migrated.provenance().rev_overlap_seq(), b"TCGT");
        assert_eq!(migrated.provenance().rev_overlap_qual(), b"IIII");
    }

    #[test]
    fn test_mergeview_lengths_and_copy_buffers_match_fixture_regions() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: "TT".into(),
            overlap_fwd_seq: "ACGT".into(),
            overlap_rev_seq: "TCGT".into(),
            right_seq: "GG".into(),
            left_qual: "II".into(),
            overlap_fwd_qual: "JKLM".into(),
            overlap_rev_qual: "WXYZ".into(),
            right_qual: "PQ".into(),
        });

        let view = validated
            .merge_view()
            .expect("merge view should construct from validated overlap");

        assert_eq!(view.left_overhang_len(), 2);
        assert_eq!(view.overlap_len(), 4);
        assert_eq!(view.right_overhang_len(), 2);
        assert_eq!(view.merged_len(), 8);

        let mut left_seq = vec![0u8; 2];
        let mut left_qual = vec![0u8; 2];
        let mut fwd_overlap_seq = vec![0u8; 4];
        let mut fwd_overlap_qual = vec![0u8; 4];
        let mut rev_overlap_seq = vec![0u8; 4];
        let mut rev_overlap_qual = vec![0u8; 4];
        let mut right_seq = vec![0u8; 2];
        let mut right_qual = vec![0u8; 2];

        view.copy_left_overhang_seq_into(&mut left_seq);
        view.copy_left_overhang_qual_into(&mut left_qual);
        view.copy_fwd_overlap_seq_into(&mut fwd_overlap_seq);
        view.copy_fwd_overlap_qual_into(&mut fwd_overlap_qual);
        view.copy_rev_overlap_seq_rc_into(&mut rev_overlap_seq);
        view.copy_rev_overlap_qual_rc_into(&mut rev_overlap_qual);
        view.copy_right_overhang_seq_rc_into(&mut right_seq);
        view.copy_right_overhang_qual_rc_into(&mut right_qual);

        assert_eq!(left_seq, b"TT");
        assert_eq!(left_qual, b"II");
        assert_eq!(fwd_overlap_seq, b"ACGT");
        assert_eq!(fwd_overlap_qual, b"JKLM");
        assert_eq!(rev_overlap_seq, b"TCGT");
        assert_eq!(rev_overlap_qual, b"WXYZ");
        assert_eq!(right_seq, b"GG");
        assert_eq!(right_qual, b"PQ");
    }

    #[test]
    fn test_copy_into_truncates_and_reports_written_len() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: "TT".into(),
            overlap_fwd_seq: "ACGT".into(),
            overlap_rev_seq: "TCGT".into(),
            right_seq: "GG".into(),
            left_qual: "II".into(),
            overlap_fwd_qual: "JKLM".into(),
            overlap_rev_qual: "WXYZ".into(),
            right_qual: "PQ".into(),
        });
        let view = validated
            .merge_view()
            .expect("merge view should construct from validated overlap");

        let mut small = [b'_'; 2];
        let written = view.copy_fwd_overlap_seq_into(&mut small);

        assert_eq!(written, 2);
        assert_eq!(&small, b"AC");
    }

    #[test]
    fn test_copy_into_oversized_buffer_writes_prefix_only() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: "TT".into(),
            overlap_fwd_seq: "ACGT".into(),
            overlap_rev_seq: "TCGT".into(),
            right_seq: "GG".into(),
            left_qual: "II".into(),
            overlap_fwd_qual: "JKLM".into(),
            overlap_rev_qual: "WXYZ".into(),
            right_qual: "PQ".into(),
        });
        let view = validated
            .merge_view()
            .expect("merge view should construct from validated overlap");

        let mut oversized = [b'_'; 6];
        let written = view.copy_right_overhang_seq_rc_into(&mut oversized);

        assert_eq!(written, 2);
        assert_eq!(&oversized[..2], b"GG");
        assert_eq!(&oversized[2..], b"____");
    }

    #[test]
    fn test_from_pair_bounds_rejects_invalid_bound_order() {
        let mates = ReadPair::from(
            SequenceRead::new("read-bounds-order", "ACGTACGT", "IIIIIIII"),
            SequenceRead::new("read-bounds-order", "ACGTACGT", "IIIIIIII"),
        )
        .expect("fixture reads should form a pair");

        let result = MergeView::from_pair_bounds(&mates, 3, 5, 3, 0, 2);
        assert!(matches!(
            result,
            Err(Error::MergeError(
                MergeError::OverlapWindowLengthMismatch { .. }
            ))
        ));
    }

    #[test]
    fn test_from_pair_bounds_rejects_forward_reverse_length_mismatch() {
        let mates = ReadPair::from(
            SequenceRead::new("read-bounds-len", "ACGTACGT", "IIIIIIII"),
            SequenceRead::new("read-bounds-len", "ACGTACGT", "IIIIIIII"),
        )
        .expect("fixture reads should form a pair");

        let result = MergeView::from_pair_bounds(&mates, 4, 2, 4, 0, 2);
        assert!(matches!(
            result,
            Err(Error::MergeError(
                MergeError::OverlapWindowLengthMismatch { .. }
            ))
        ));
    }

    proptest! {
        #[test]
        fn proptest_merge_matches_oracle_for_constructed_overlap(
            fixture in merge_fixture_strategy(),
        ) {
            let validated = build_validated_overlap_from_fixture(&fixture);

            let merged = merge_from(&validated)
                .expect("constructed overlap fixture should merge via generic entrypoint");
            let (expected_seq, expected_qual) = oracle_merge(&fixture);

            prop_assert_eq!(merged.sequence(), expected_seq.as_slice());
            prop_assert_eq!(merged.qualities(), expected_qual.as_slice());
            prop_assert_eq!(merged.sequence().len(), merged.qualities().len());
            prop_assert_eq!(
                merged.provenance().fwd_overlap_seq(),
                fixture.overlap_fwd_seq.as_bytes()
            );
            prop_assert_eq!(
                merged.provenance().fwd_overlap_qual(),
                fixture.overlap_fwd_qual.as_bytes()
            );
            prop_assert_eq!(
                merged.provenance().rev_overlap_seq(),
                fixture.overlap_rev_seq.as_bytes()
            );
            prop_assert_eq!(
                merged.provenance().rev_overlap_qual(),
                fixture.overlap_rev_qual.as_bytes()
            );
        }
    }
}
