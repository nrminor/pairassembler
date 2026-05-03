use std::{cmp::Ordering, ops::Range};

use crate::{
    Result,
    assembler::HasPairOverlap,
    errors::MergeError::{
        EmptyOverlapWindow, EqualQualityBaseDisagreement, MergeSequenceQualityLengthMismatch,
        MergedLengthMismatch, OverlapWindowLengthMismatch, ProvenanceLengthMismatch,
    },
    overlap::{HasOrientedPairSlices, OverlapBounds},
    prelude::utils::encode_fastq_quality_scores_in_place,
    read::OwnedSequenceRead,
};

/// Parameters controlling deterministic overlap merging.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MergeParams {
    tie_policy: MergeTiePolicy,
}

/// Policy for equal-quality, disagreeing base calls inside an already-selected overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MergeTiePolicy {
    /// Choose the forward mate's base.
    #[default]
    PreferForward,
    /// Choose the reverse mate's base in reverse-complement orientation.
    PreferReverse,
    /// Emit `N` rather than pretending either equal-quality disagreement is better supported.
    EmitAmbiguous,
    /// Fail merge when equal-quality base calls disagree.
    RejectDisagreement,
    /// Prefer the base farther from its source read end, falling back to forward on distance ties.
    PreferInteriorBase,
}

/// Applies overlap merge policy to pair-overlap slices.
pub(crate) struct OverlapMerger {
    params: MergeParams,
}

impl MergeParams {
    #[must_use]
    pub fn with_tie_policy(self, tie_policy: MergeTiePolicy) -> Self {
        Self { tie_policy }
    }

    #[must_use]
    pub(crate) fn tie_policy(self) -> MergeTiePolicy {
        self.tie_policy
    }
}

/// Score-space consensus produced by the merge kernel.
///
/// This is the minimal merged payload needed by staged contexts. It deliberately does not own
/// overlap provenance; contexts that still need correction should retain the overlap slices they
/// already carried into merge.
///
/// `MergedConsensus` is for staged contexts that still retain `PairOverlap` separately.
/// `MergedRead` is for detached merged reads and therefore owns overlap provenance.
#[derive(Debug, Clone)]
pub(crate) struct MergedConsensus {
    pub(crate) id: String,
    pub(crate) sequence: Vec<u8>,
    pub(crate) quality_scores: Vec<u8>,
    pub(crate) left_overhang_len: usize,
}

/// Deterministic consensus read produced by the merge stage.
///
/// This is the canonical merge artifact for downstream processing.
#[derive(Debug, Clone)]
pub struct MergedRead {
    pub(crate) id: String,
    pub(crate) consensus_seq: Vec<u8>,
    pub(crate) consensus_quality_scores: Vec<u8>,
    pub(crate) left_overhang_len: usize,
    pub(crate) provenance: MergeProvenance,
}

/// Overlap-only windows retained from merging.
///
/// This payload is intentionally narrow: it records only overlap windows needed
/// for downstream correction/diagnostics and does not represent full-read
/// provenance.
#[derive(Debug, Clone)]
pub(crate) struct MergeProvenance {
    overlap_len: usize,
    fwd_overlap_seq: Vec<u8>,
    fwd_overlap_quality_scores: Vec<u8>,
    rev_overlap_seq: Vec<u8>,
    rev_overlap_quality_scores: Vec<u8>,
}

#[derive(Debug, Default)]
pub(crate) struct MergeProvenanceBuilder {
    fwd_overlap_seq: Option<Vec<u8>>,
    fwd_overlap_quality_scores: Option<Vec<u8>>,
    rev_overlap_seq: Option<Vec<u8>>,
    rev_overlap_quality_scores: Option<Vec<u8>>,
}

/// Borrowed normalized merge input projected from a pair-plus-overlap carrier.
///
/// `MergeView` is an internal normalization boundary rather than a domain object:
/// it packages the exact slices and coordinate mappings needed by the cheap merge
/// kernel without owning sequence data.
#[derive(Debug, Clone, Copy)]
struct MergeView<'a> {
    id: &'a str,
    fwd_overlap_start: usize,
    fwd_read_len: usize,
    rev_overlap_start: usize,
    rev_read_len: usize,
    left_overhang_seq: &'a [u8],
    left_overhang_qual: &'a [u8],
    fwd_overlap_seq: &'a [u8],
    fwd_overlap_qual: &'a [u8],
    rev_overlap_seq_rc: &'a [u8],
    rev_overlap_qual_rc: &'a [u8],
    right_overhang_seq_rc: &'a [u8],
    right_overhang_qual_rc: &'a [u8],
}

#[derive(Debug, Clone)]
struct CheckedOverlapRanges {
    fwd: Range<usize>,
    rev_rc: Range<usize>,
}

impl<'a> MergeView<'a> {
    fn from_pair_overlap<T>(input: &'a T) -> Result<Self>
    where
        T: HasPairOverlap + ?Sized,
    {
        input.validate_overlap_bounds()?;
        let slices = input.pair_slices()?;
        let bounds = input.overlap_bounds()?;

        let (fwd_seq, rev_seq_rc) = slices.sequences();
        let (fwd_qual, rev_qual_rc) = slices.quality_score_bytes();

        let ranges = CheckedOverlapRanges::from_bounds(bounds, fwd_seq.len(), rev_seq_rc.len())?;

        ensure_seq_qual_lengths("forward", fwd_seq, fwd_qual)?;
        ensure_seq_qual_lengths("reverse_rc", rev_seq_rc, rev_qual_rc)?;

        Ok(Self {
            id: slices.pair_id(),
            fwd_overlap_start: ranges.fwd.start,
            fwd_read_len: fwd_seq.len(),
            rev_overlap_start: ranges.rev_rc.start,
            rev_read_len: rev_seq_rc.len(),
            left_overhang_seq: &fwd_seq[..ranges.fwd.start],
            left_overhang_qual: &fwd_qual[..ranges.fwd.start],
            fwd_overlap_seq: &fwd_seq[ranges.fwd.clone()],
            fwd_overlap_qual: &fwd_qual[ranges.fwd],
            rev_overlap_seq_rc: &rev_seq_rc[ranges.rev_rc.clone()],
            rev_overlap_qual_rc: &rev_qual_rc[ranges.rev_rc.clone()],
            right_overhang_seq_rc: &rev_seq_rc[ranges.rev_rc.end..],
            right_overhang_qual_rc: &rev_qual_rc[ranges.rev_rc.end..],
        })
    }

    #[inline]
    fn id(&self) -> &str {
        self.id
    }

    #[inline]
    fn left_overhang_len(&self) -> usize {
        self.left_overhang_seq.len()
    }

    #[inline]
    fn overlap_len(&self) -> usize {
        self.fwd_overlap_seq.len()
    }

    #[inline]
    fn right_overhang_len(&self) -> usize {
        self.right_overhang_seq_rc.len()
    }

    #[inline]
    fn merged_len(&self) -> usize {
        self.left_overhang_len() + self.overlap_len() + self.right_overhang_len()
    }

    #[inline]
    fn left_overhang_seq(&self) -> &[u8] {
        self.left_overhang_seq
    }

    #[inline]
    fn left_overhang_qual(&self) -> &[u8] {
        self.left_overhang_qual
    }

    #[inline]
    fn fwd_overlap_seq(&self) -> &[u8] {
        self.fwd_overlap_seq
    }

    #[inline]
    fn fwd_overlap_qual(&self) -> &[u8] {
        self.fwd_overlap_qual
    }

    #[inline]
    fn overlap_pair_at(&self, i: usize) -> ((u8, u8), (u8, u8)) {
        debug_assert!(i < self.overlap_len());
        let fwd = (self.fwd_overlap_seq[i], self.fwd_overlap_qual()[i]);
        let rev = (self.rev_overlap_seq_rc[i], self.rev_overlap_qual_rc[i]);
        (fwd, rev)
    }

    #[inline]
    fn overlap_interior_distances_at(&self, i: usize) -> (usize, usize) {
        debug_assert!(i < self.overlap_len());
        let fwd_pos = self.fwd_overlap_start + i;
        let rev_pos = self.rev_overlap_start + i;
        (
            distance_to_nearest_read_end(fwd_pos, self.fwd_read_len),
            distance_to_nearest_read_end(rev_pos, self.rev_read_len),
        )
    }

    #[inline]
    fn right_overhang_pair_at(&self, i: usize) -> (u8, u8) {
        (
            self.right_overhang_seq_rc[i],
            self.right_overhang_qual_rc[i],
        )
    }
}

impl CheckedOverlapRanges {
    fn from_bounds(bounds: OverlapBounds, fwd_len: usize, rev_rc_len: usize) -> Result<Self> {
        let overlap_len = bounds.overlap_len();
        if overlap_len == 0 {
            return Err(EmptyOverlapWindow.into());
        }

        let fwd = bounds.forward_range();
        let rev_rc = bounds.reverse_range();

        if fwd.start >= fwd.end
            || rev_rc.start >= rev_rc.end
            || fwd.end > fwd_len
            || rev_rc.end > rev_rc_len
        {
            return Err(OverlapWindowLengthMismatch {
                fwd_len: fwd.end.saturating_sub(fwd.start),
                rev_len: rev_rc.end.saturating_sub(rev_rc.start),
            }
            .into());
        }

        if fwd.len() != overlap_len || rev_rc.len() != overlap_len {
            return Err(OverlapWindowLengthMismatch {
                fwd_len: fwd.len(),
                rev_len: rev_rc.len(),
            }
            .into());
        }

        Ok(Self { fwd, rev_rc })
    }
}

impl MergedRead {
    /// Construct a merged read from checked consensus and retained merge provenance.
    ///
    /// # Errors
    ///
    /// Returns an error if consensus sequence and quality lengths differ, or if the retained
    /// overlap window cannot fit within the consensus layout described by `left_overhang_len`.
    pub(crate) fn from_consensus_and_provenance(
        consensus: MergedConsensus,
        provenance: MergeProvenance,
    ) -> Result<Self> {
        let MergedConsensus {
            id,
            sequence: consensus_seq,
            quality_scores: consensus_quality_scores,
            left_overhang_len,
        } = consensus;

        ensure_seq_qual_lengths("consensus", &consensus_seq, &consensus_quality_scores)?;
        if left_overhang_len + provenance.overlap_len() > consensus_seq.len() {
            return Err(MergedLengthMismatch {
                expected: left_overhang_len + provenance.overlap_len(),
                actual: consensus_seq.len(),
            }
            .into());
        }

        Ok(Self {
            id,
            consensus_seq,
            consensus_quality_scores,
            left_overhang_len,
            provenance,
        })
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

    /// Borrow merged consensus quality score bytes.
    #[must_use]
    pub fn quality_score_bytes(&self) -> &[u8] {
        self.consensus_quality_scores.as_slice()
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

    /// Consume and return owned merged quality score bytes.
    #[must_use]
    pub fn quality_score_bytes_owned(self) -> Vec<u8> {
        self.consensus_quality_scores
    }

    /// Return ASCII-encoded FASTQ quality bytes for the merged consensus.
    #[must_use]
    pub fn to_quality_ascii_bytes(&self) -> Vec<u8> {
        let mut quality_ascii = self.consensus_quality_scores.clone();
        encode_fastq_quality_scores_in_place(&mut quality_ascii);
        quality_ascii
    }

    /// Borrow overlap windows retained from merge.
    #[cfg(test)]
    #[must_use]
    fn provenance(&self) -> &MergeProvenance {
        &self.provenance
    }
}

impl MergedConsensus {
    /// Construct checked merged consensus parts.
    ///
    /// # Errors
    ///
    /// Returns an error if consensus sequence and quality score lengths differ.
    pub(crate) fn try_new(
        id: String,
        sequence: Vec<u8>,
        quality_scores: Vec<u8>,
        left_overhang_len: usize,
    ) -> Result<Self> {
        ensure_seq_qual_lengths("consensus", &sequence, &quality_scores)?;

        Ok(Self {
            id,
            sequence,
            quality_scores,
            left_overhang_len,
        })
    }

    /// Borrow the merged consensus identifier.
    #[must_use]
    pub(crate) fn id(&self) -> &str {
        &self.id
    }

    /// Borrow merged consensus sequence bytes.
    #[must_use]
    pub(crate) fn sequence(&self) -> &[u8] {
        self.sequence.as_slice()
    }

    /// Borrow merged consensus quality score bytes.
    #[must_use]
    pub(crate) fn quality_score_bytes(&self) -> &[u8] {
        self.quality_scores.as_slice()
    }

    /// Return the consensus left-overhang length before the overlap window.
    #[must_use]
    pub(crate) fn left_overhang_len(&self) -> usize {
        self.left_overhang_len
    }
}

impl OverlapMerger {
    pub(crate) fn new(params: MergeParams) -> Self {
        Self { params }
    }

    pub(crate) fn merge_pair_overlap<T>(&self, input: &T) -> Result<MergedRead>
    where
        T: HasPairOverlap + ?Sized,
    {
        self.merge_view(MergeView::from_pair_overlap(input)?)
    }

    pub(crate) fn merge_consensus<T>(&self, input: &T) -> Result<MergedConsensus>
    where
        T: HasPairOverlap + ?Sized,
    {
        let view = MergeView::from_pair_overlap(input)?;
        self.merge_consensus_view(&view)
    }

    fn merge_view(&self, view: MergeView<'_>) -> Result<MergedRead> {
        let overlap_len = view.overlap_len();
        debug_assert_eq!(overlap_len, view.fwd_overlap_seq().len());

        let consensus = self.merge_consensus_view(&view)?;

        let provenance = MergeProvenance::builder()
            .forward_overlap(
                view.fwd_overlap_seq().to_vec(),
                view.fwd_overlap_qual().to_vec(),
            )
            .reverse_overlap_rc(
                view.rev_overlap_seq_rc.to_vec(),
                view.rev_overlap_qual_rc.to_vec(),
            )
            .build()?;

        MergedRead::from_consensus_and_provenance(consensus, provenance)
    }

    fn merge_consensus_view(&self, view: &MergeView<'_>) -> Result<MergedConsensus> {
        let overlap_len = view.overlap_len();
        debug_assert_eq!(overlap_len, view.fwd_overlap_seq().len());

        let expected_len = view.merged_len();
        let mut full_seq = Vec::with_capacity(expected_len);
        let mut full_qual = Vec::with_capacity(expected_len);

        full_seq.extend_from_slice(view.left_overhang_seq());
        full_qual.extend_from_slice(view.left_overhang_qual());

        for i in 0..overlap_len {
            let (base, quality) = self.merge_overlap_column(view, i)?;
            full_seq.push(base);
            full_qual.push(quality);
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

        MergedConsensus::try_new(
            view.id().to_owned(),
            full_seq,
            full_qual,
            view.left_overhang_len(),
        )
    }

    fn merge_overlap_column(&self, view: &MergeView<'_>, offset: usize) -> Result<(u8, u8)> {
        let ((fwd_base, fwd_quality), (rev_base, rev_quality)) = view.overlap_pair_at(offset);

        match fwd_quality.cmp(&rev_quality) {
            Ordering::Greater => Ok((fwd_base, fwd_quality)),
            Ordering::Less => Ok((rev_base, rev_quality)),
            Ordering::Equal => {
                self.merge_equal_quality_column(view, offset, fwd_base, rev_base, fwd_quality)
            },
        }
    }

    fn merge_equal_quality_column(
        &self,
        view: &MergeView<'_>,
        offset: usize,
        fwd_base: u8,
        rev_base: u8,
        quality: u8,
    ) -> Result<(u8, u8)> {
        if fwd_base == rev_base {
            return Ok((fwd_base, quality));
        }

        match self.params.tie_policy() {
            MergeTiePolicy::PreferForward => Ok((fwd_base, quality)),
            MergeTiePolicy::PreferReverse => Ok((rev_base, quality)),
            MergeTiePolicy::EmitAmbiguous => Ok((b'N', quality)),
            MergeTiePolicy::RejectDisagreement => Err(EqualQualityBaseDisagreement {
                offset,
                fwd_base,
                rev_base,
                quality,
            }
            .into()),
            MergeTiePolicy::PreferInteriorBase => {
                let (fwd_distance, rev_distance) = view.overlap_interior_distances_at(offset);
                if fwd_distance >= rev_distance {
                    Ok((fwd_base, quality))
                } else {
                    Ok((rev_base, quality))
                }
            },
        }
    }
}

impl MergeProvenance {
    pub(crate) fn builder() -> MergeProvenanceBuilder {
        MergeProvenanceBuilder::default()
    }

    fn from_builder(
        fwd_overlap_seq: Vec<u8>,
        fwd_overlap_quality_scores: Vec<u8>,
        rev_overlap_seq: Vec<u8>,
        rev_overlap_quality_scores: Vec<u8>,
    ) -> Result<Self> {
        let overlap_len = fwd_overlap_seq.len();
        ensure_seq_qual_lengths("overlap_fwd", &fwd_overlap_seq, &fwd_overlap_quality_scores)?;
        ensure_seq_qual_lengths("overlap_rev", &rev_overlap_seq, &rev_overlap_quality_scores)?;

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
            fwd_overlap_quality_scores,
            rev_overlap_seq,
            rev_overlap_quality_scores,
        })
    }

    /// Return overlap length used by merge.
    #[must_use]
    pub(crate) fn overlap_len(&self) -> usize {
        self.overlap_len
    }

    /// Borrow forward overlap sequence bytes.
    #[must_use]
    pub(crate) fn fwd_overlap_seq(&self) -> &[u8] {
        self.fwd_overlap_seq.as_slice()
    }

    /// Borrow forward overlap quality score bytes.
    #[must_use]
    pub(crate) fn fwd_overlap_quality_score_bytes(&self) -> &[u8] {
        self.fwd_overlap_quality_scores.as_slice()
    }

    /// Borrow reverse overlap sequence bytes in reverse-complement orientation.
    #[must_use]
    pub(crate) fn rev_overlap_seq(&self) -> &[u8] {
        self.rev_overlap_seq.as_slice()
    }

    /// Borrow reverse overlap quality score bytes in reverse-complement orientation.
    #[must_use]
    pub(crate) fn rev_overlap_quality_score_bytes(&self) -> &[u8] {
        self.rev_overlap_quality_scores.as_slice()
    }
}

impl MergeProvenanceBuilder {
    pub(crate) fn forward_overlap(mut self, seq: Vec<u8>, quality_scores: Vec<u8>) -> Self {
        self.fwd_overlap_seq = Some(seq);
        self.fwd_overlap_quality_scores = Some(quality_scores);
        self
    }

    pub(crate) fn reverse_overlap_rc(mut self, seq: Vec<u8>, quality_scores: Vec<u8>) -> Self {
        self.rev_overlap_seq = Some(seq);
        self.rev_overlap_quality_scores = Some(quality_scores);
        self
    }

    pub(crate) fn build(self) -> Result<MergeProvenance> {
        let fwd_overlap_seq = Self::required(self.fwd_overlap_seq, "forward overlap sequence")?;
        let fwd_overlap_quality_scores = Self::required(
            self.fwd_overlap_quality_scores,
            "forward overlap quality scores",
        )?;
        let rev_overlap_seq = Self::required(self.rev_overlap_seq, "reverse overlap sequence")?;
        let rev_overlap_quality_scores = Self::required(
            self.rev_overlap_quality_scores,
            "reverse overlap quality scores",
        )?;

        MergeProvenance::from_builder(
            fwd_overlap_seq,
            fwd_overlap_quality_scores,
            rev_overlap_seq,
            rev_overlap_quality_scores,
        )
    }

    fn required<T>(value: Option<T>, name: &'static str) -> Result<T> {
        value.ok_or_else(|| anyhow::anyhow!("missing {name} for merge provenance").into())
    }
}

impl TryFrom<MergedRead> for OwnedSequenceRead {
    type Error = crate::Error;

    fn try_from(read: MergedRead) -> Result<Self> {
        let mut quality_ascii = read.consensus_quality_scores;
        encode_fastq_quality_scores_in_place(&mut quality_ascii);

        Self::try_from_ascii_bytes(read.id, read.consensus_seq, quality_ascii)
    }
}

impl TryFrom<MergedConsensus> for OwnedSequenceRead {
    type Error = crate::Error;

    fn try_from(consensus: MergedConsensus) -> Result<Self> {
        let mut quality_ascii = consensus.quality_scores;
        encode_fastq_quality_scores_in_place(&mut quality_ascii);

        Self::try_from_ascii_bytes(consensus.id, consensus.sequence, quality_ascii)
    }
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

#[inline]
fn distance_to_nearest_read_end(position: usize, read_len: usize) -> usize {
    debug_assert!(position < read_len);
    position.min(read_len - 1 - position)
}

#[cfg(test)]
mod tests {
    use super::{CheckedOverlapRanges, MergeParams, MergeTiePolicy, MergeView, OverlapMerger};
    use crate::{
        Error, PairOverlap, Result,
        errors::MergeError,
        overlap::{OrientedPairSlices, OverlapBounds},
        prelude::utils::decode_fastq_quality_scores,
        read::{ReadPair, SequenceRead},
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

    fn overlap_from_mates<'a>(
        mates: &'a ReadPair<'a>,
        overlap_len: usize,
        fwd_start_offset: usize,
        rev_start_offset: usize,
    ) -> PairOverlap<'a> {
        let slices = mates.to_oriented_slices();

        PairOverlap::from_oriented_slices(
            slices,
            OverlapBounds::new(overlap_len, fwd_start_offset, rev_start_offset),
        )
        .expect("test overlap should satisfy overlap invariants")
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
        let overlap = overlap_from_mates(mates_ref, overlap_len, fwd_start, 0);

        let metrics = ValidationMetrics::new(overlap_len, overlap_len, 0, 0.0);

        ValidatedOverlap::new_unchecked(overlap, metrics)
    }

    fn merge_with_default_params(validated: &ValidatedOverlap<'_>) -> Result<super::MergedRead> {
        OverlapMerger::new(MergeParams::default()).merge_pair_overlap(validated)
    }

    fn merge_with_tie_policy(
        validated: &ValidatedOverlap<'_>,
        tie_policy: MergeTiePolicy,
    ) -> Result<super::MergedRead> {
        OverlapMerger::new(MergeParams::default().with_tie_policy(tie_policy))
            .merge_pair_overlap(validated)
    }

    fn oracle_merge(fixture: &MergeFixture) -> (Vec<u8>, Vec<u8>) {
        let left_seq = fixture.left_seq.as_str();
        let overlap_fwd_seq = fixture.overlap_fwd_seq.as_str();
        let overlap_rev_seq = fixture.overlap_rev_seq.as_str();
        let right_seq = fixture.right_seq.as_str();
        let left_qual = decode_fastq_quality_scores(fixture.left_qual.as_bytes());
        let overlap_fwd_qual = decode_fastq_quality_scores(fixture.overlap_fwd_qual.as_bytes());
        let overlap_rev_qual = decode_fastq_quality_scores(fixture.overlap_rev_qual.as_bytes());
        let right_qual = decode_fastq_quality_scores(fixture.right_qual.as_bytes());

        let mut overlap_seq = Vec::with_capacity(overlap_fwd_seq.len());
        let mut overlap_qual = Vec::with_capacity(overlap_fwd_qual.len());

        for (((fwd_base, fwd_q), rev_base), rev_q) in overlap_fwd_seq
            .bytes()
            .zip(overlap_fwd_qual.iter().copied())
            .zip(overlap_rev_seq.bytes())
            .zip(overlap_rev_qual.iter().copied())
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
        full_qual.extend_from_slice(&left_qual);
        full_qual.extend_from_slice(&overlap_qual);
        full_qual.extend_from_slice(&right_qual);

        (full_seq, full_qual)
    }

    #[test]
    fn test_merge_perfect_full_overlap_roundtrip() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = overlap_from_mates(&mates, 5, 4, 0);
        let metrics = ValidationMetrics::new(5, 5, 0, 0.0);
        let validated = ValidatedOverlap::new_unchecked(overlap, metrics);

        let merged = merge_with_default_params(&validated)
            .expect("overlap merger should merge validated overlap without bounds errors");

        assert_eq!(merged.id(), "read1");
        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.quality_score_bytes(), &[40; 9]);
        assert_eq!(merged.to_quality_ascii_bytes(), b"IIIIIIIII");
        assert_eq!(merged.sequence().len(), merged.quality_score_bytes().len());
    }

    #[test]
    fn test_merge_with_left_overhang_preserves_prefix() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = overlap_from_mates(&mates, 5, 4, 0);
        let metrics = ValidationMetrics::new(5, 5, 0, 0.0);
        let validated = ValidatedOverlap::new_unchecked(overlap, metrics);

        let merged = merge_with_default_params(&validated)
            .expect("overlap merger should merge validated overlap without bounds errors");

        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.sequence().len(), merged.quality_score_bytes().len());
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

        let merged = merge_with_default_params(&validated)
            .expect("overlap merger should preserve right overhang semantics");

        assert_eq!(merged.sequence(), b"TTACGTGG");
        assert_eq!(merged.quality_score_bytes(), &[40; 8]);
        assert_eq!(merged.to_quality_ascii_bytes(), b"IIIIIIII");
        assert_eq!(merged.sequence().len(), merged.quality_score_bytes().len());
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

        let merged = merge_with_default_params(&validated)
            .expect("overlap merger should remain deterministic on equal-quality overlap");

        assert_eq!(merged.sequence(), b"AAAA");
        assert_eq!(merged.quality_score_bytes(), &[40; 4]);
    }

    #[test]
    fn test_merge_tie_on_quality_can_prefer_reverse_base() {
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

        let merged = merge_with_tie_policy(&validated, MergeTiePolicy::PreferReverse)
            .expect("reverse tie policy should merge equal-quality disagreements");

        assert_eq!(merged.sequence(), b"TTTT");
        assert_eq!(merged.quality_score_bytes(), &[40; 4]);
    }

    #[test]
    fn test_merge_tie_on_quality_can_emit_ambiguous_base() {
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

        let merged = merge_with_tie_policy(&validated, MergeTiePolicy::EmitAmbiguous)
            .expect("ambiguous tie policy should merge equal-quality disagreements");

        assert_eq!(merged.sequence(), b"NNNN");
        assert_eq!(merged.quality_score_bytes(), &[40; 4]);
    }

    #[test]
    fn test_merge_tie_on_quality_can_reject_disagreement() {
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

        let err = merge_with_tie_policy(&validated, MergeTiePolicy::RejectDisagreement)
            .expect_err("reject tie policy should fail equal-quality disagreements");

        assert!(matches!(
            err,
            Error::MergeError(MergeError::EqualQualityBaseDisagreement {
                offset: 0,
                fwd_base: b'A',
                rev_base: b'T',
                quality: 40,
            })
        ));
    }

    #[test]
    fn test_merge_tie_on_quality_can_prefer_interior_base() {
        let validated = build_validated_overlap_from_fixture(&MergeFixture {
            left_seq: String::new(),
            overlap_fwd_seq: "AAAA".into(),
            overlap_rev_seq: "TTTT".into(),
            right_seq: "GGGG".into(),
            left_qual: String::new(),
            overlap_fwd_qual: "IIII".into(),
            overlap_rev_qual: "IIII".into(),
            right_qual: "IIII".into(),
        });

        let merged = merge_with_tie_policy(&validated, MergeTiePolicy::PreferInteriorBase)
            .expect("interior-base tie policy should merge equal-quality disagreements");

        assert_eq!(merged.sequence(), b"AATTGGGG");
        assert_eq!(merged.quality_score_bytes(), &[40; 8]);
    }

    #[test]
    fn test_overlap_merger_populates_provenance_overlap_windows() {
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

        let migrated = merge_with_default_params(&validated)
            .expect("overlap merger should succeed on validated overlap fixture");

        assert_eq!(migrated.provenance().overlap_len(), 4);
        assert_eq!(migrated.provenance().fwd_overlap_seq(), b"ACGT");
        assert_eq!(
            migrated.provenance().fwd_overlap_quality_score_bytes(),
            &[40; 4]
        );
        assert_eq!(migrated.provenance().rev_overlap_seq(), b"TCGT");
        assert_eq!(
            migrated.provenance().rev_overlap_quality_score_bytes(),
            &[40; 4]
        );
    }

    #[test]
    fn test_mergeview_lengths_and_regions_match_fixture() {
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

        let view = MergeView::from_pair_overlap(&validated)
            .expect("merge view should construct from validated overlap");

        assert_eq!(view.left_overhang_len(), 2);
        assert_eq!(view.overlap_len(), 4);
        assert_eq!(view.right_overhang_len(), 2);
        assert_eq!(view.merged_len(), 8);

        let rev_overlap: Vec<_> = (0..view.overlap_len())
            .map(|offset| view.overlap_pair_at(offset).1.0)
            .collect();
        let rev_overlap_qual: Vec<_> = (0..view.overlap_len())
            .map(|offset| view.overlap_pair_at(offset).1.1)
            .collect();
        let right: Vec<_> = (0..view.right_overhang_len())
            .map(|offset| view.right_overhang_pair_at(offset).0)
            .collect();
        let right_qual: Vec<_> = (0..view.right_overhang_len())
            .map(|offset| view.right_overhang_pair_at(offset).1)
            .collect();

        assert_eq!(view.left_overhang_seq(), b"TT");
        assert_eq!(view.left_overhang_qual(), [40, 40]);
        assert_eq!(view.fwd_overlap_seq(), b"ACGT");
        assert_eq!(view.fwd_overlap_qual(), [41, 42, 43, 44]);
        assert_eq!(rev_overlap, b"TCGT");
        assert_eq!(rev_overlap_qual, [54, 55, 56, 57]);
        assert_eq!(right, b"GG");
        assert_eq!(right_qual, [47, 48]);
    }

    #[test]
    fn test_checked_overlap_ranges_rejects_forward_out_of_bounds() {
        let result = CheckedOverlapRanges::from_bounds(OverlapBounds::new(4, 5, 0), 8, 8);

        assert!(matches!(
            result,
            Err(Error::MergeError(
                MergeError::OverlapWindowLengthMismatch { .. }
            ))
        ));
    }

    #[test]
    fn test_checked_overlap_ranges_rejects_reverse_out_of_bounds() {
        let result = CheckedOverlapRanges::from_bounds(OverlapBounds::new(4, 2, 0), 8, 3);

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

            let merged = merge_with_default_params(&validated)
                .expect("constructed overlap fixture should merge via overlap merger");
            let (expected_seq, expected_qual) = oracle_merge(&fixture);

            prop_assert_eq!(merged.sequence(), expected_seq.as_slice());
            prop_assert_eq!(merged.quality_score_bytes(), expected_qual.as_slice());
            prop_assert_eq!(merged.sequence().len(), merged.quality_score_bytes().len());
            prop_assert_eq!(
                merged.provenance().fwd_overlap_seq(),
                fixture.overlap_fwd_seq.as_bytes()
            );
            let expected_fwd_overlap_quality_scores =
                decode_fastq_quality_scores(fixture.overlap_fwd_qual.as_bytes());
            prop_assert_eq!(
                merged.provenance().fwd_overlap_quality_score_bytes(),
                expected_fwd_overlap_quality_scores.as_ref()
            );
            prop_assert_eq!(
                merged.provenance().rev_overlap_seq(),
                fixture.overlap_rev_seq.as_bytes()
            );
            let expected_rev_overlap_quality_scores =
                decode_fastq_quality_scores(fixture.overlap_rev_qual.as_bytes());
            prop_assert_eq!(
                merged.provenance().rev_overlap_quality_score_bytes(),
                expected_rev_overlap_quality_scores.as_ref()
            );
        }
    }
}
