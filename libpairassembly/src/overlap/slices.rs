use crate::{Result, errors::OverlapError, prelude::utils::fastq_ascii_to_phred, read::ReadPair};

use super::{OverlapBounds, private};
#[cfg(test)]
use super::{OverlapFinder, PairOverlap};

#[allow(clippy::elidable_lifetime_names)]
impl<'a> ReadPair<'a> {
    /// Discover the best overlap between this read pair.
    ///
    /// Returns `Ok(None)` when no overlap candidate satisfies the configured thresholds.
    ///
    /// # Errors
    ///
    /// Returns an error when overlap candidate reconciliation fails (for example, tie rejection),
    /// or when computed overlap coordinates are inconsistent with read bounds.
    #[cfg(test)]
    pub(crate) fn overlap<'scratch>(
        &self,
        params: &super::OverlapParams,
        scratch: &'scratch mut AssemblyScratch,
    ) -> Result<Option<PairOverlap<'a, 'scratch>>> {
        OverlapFinder::new(params).find(*self, scratch)
    }

    #[cfg(test)]
    pub(crate) fn to_oriented_slices(
        self,
        scratch: &mut AssemblyScratch,
    ) -> OrientedPairSlices<'a, '_> {
        scratch.orient(self)
    }
}

#[derive(Debug, Default)]
pub(crate) struct AssemblyScratch {
    fwd_quality_score_bytes: Vec<u8>,
    rev_seq_rc: Vec<u8>,
    rev_quality_score_bytes_rc: Vec<u8>,
}

impl AssemblyScratch {
    pub(crate) fn orient<'pair, 'scratch>(
        &'scratch mut self,
        pair: ReadPair<'pair>,
    ) -> OrientedPairSlices<'pair, 'scratch> {
        self.fwd_quality_score_bytes.clear();
        self.fwd_quality_score_bytes.extend(
            pair.fwd_quality_bytes()
                .iter()
                .copied()
                .map(fastq_ascii_to_phred),
        );

        self.rev_seq_rc.clear();
        self.rev_seq_rc.extend(
            pair.rev_sequence_bytes()
                .iter()
                .rev()
                .copied()
                .map(complement_base),
        );

        self.rev_quality_score_bytes_rc.clear();
        self.rev_quality_score_bytes_rc.extend(
            pair.rev_quality_bytes()
                .iter()
                .rev()
                .copied()
                .map(fastq_ascii_to_phred),
        );

        OrientedPairSlices {
            id: pair.fwd_id(),
            fwd_seq: pair.fwd_sequence_bytes(),
            fwd_quality_score_bytes: &self.fwd_quality_score_bytes,
            rev_seq_rc: &self.rev_seq_rc,
            rev_quality_score_bytes_rc: &self.rev_quality_score_bytes_rc,
        }
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
/// Retained pair slices after converting reads into overlap orientation.
///
/// Forward sequence bytes are borrowed from the original read. Quality scores and reverse-mate
/// orientation borrow buffers from the owning assembler scratch space because FASTQ ASCII qualities
/// must be decoded and the reverse mate must be reverse-complemented before downstream validation,
/// merge, and correction.
pub struct OrientedPairSlices<'pair, 'scratch> {
    pub(crate) id: &'pair str,
    pub(crate) fwd_seq: &'pair [u8],
    pub(crate) fwd_quality_score_bytes: &'scratch [u8],
    pub(crate) rev_seq_rc: &'scratch [u8],
    pub(crate) rev_quality_score_bytes_rc: &'scratch [u8],
}

/// Exposes score-space paired slices in overlap orientation.
///
/// Implementors must preserve these invariants:
/// - forward sequence is in forward-read orientation;
/// - forward qualities are numeric quality score bytes;
/// - reverse sequence is reverse-complemented into forward-read orientation;
/// - reverse qualities are reversed to match the reverse-complemented sequence;
/// - sequence and quality lengths match within each mate.
pub(crate) trait HasOrientedPairSlices: private::Sealed {
    fn pair_id(&self) -> &str;
    fn forward_sequence(&self) -> &[u8];
    fn forward_quality_score_bytes(&self) -> &[u8];
    fn reverse_sequence_rc(&self) -> &[u8];
    fn reverse_quality_score_bytes_rc(&self) -> &[u8];

    fn forward_len(&self) -> usize {
        self.forward_sequence().len()
    }

    fn reverse_len(&self) -> usize {
        self.reverse_sequence_rc().len()
    }

    fn sequences(&self) -> (&[u8], &[u8]) {
        (self.forward_sequence(), self.reverse_sequence_rc())
    }

    fn quality_score_bytes(&self) -> (&[u8], &[u8]) {
        (
            self.forward_quality_score_bytes(),
            self.reverse_quality_score_bytes_rc(),
        )
    }

    fn validate_shape(&self) -> Result<()> {
        ensure_sequence_quality_lengths(
            "fwd_mate",
            self.forward_sequence().len(),
            self.forward_quality_score_bytes().len(),
        )?;
        ensure_sequence_quality_lengths(
            "rev_mate_rc",
            self.reverse_sequence_rc().len(),
            self.reverse_quality_score_bytes_rc().len(),
        )
    }

    fn validate_overlap_bounds(&self, bounds: OverlapBounds) -> Result<()> {
        self.validate_shape()?;

        if bounds.overlap_len() == 0 {
            return Err(OverlapError::InvalidOverlapLength {
                computed: bounds.overlap_len(),
                read1_len: self.forward_len(),
                read2_len: self.reverse_len(),
                min_required: 1,
            }
            .into());
        }

        let fwd_end = bounds.fwd_end_offset();
        let rev_end = bounds.rev_end_offset();

        let fwd_len = self.forward_len();
        if fwd_end >= fwd_len {
            return Err(OverlapError::IndexOutOfBounds {
                read: "fwd_mate",
                index: fwd_end,
                length: fwd_len,
            }
            .into());
        }

        let rev_len = self.reverse_len();
        if rev_end >= rev_len {
            return Err(OverlapError::IndexOutOfBounds {
                read: "rev_mate",
                index: rev_end,
                length: rev_len,
            }
            .into());
        }

        Ok(())
    }
}

impl private::Sealed for OrientedPairSlices<'_, '_> {}

impl HasOrientedPairSlices for OrientedPairSlices<'_, '_> {
    fn pair_id(&self) -> &str {
        self.id
    }

    fn forward_sequence(&self) -> &[u8] {
        self.fwd_seq
    }

    fn forward_quality_score_bytes(&self) -> &[u8] {
        self.fwd_quality_score_bytes
    }

    fn reverse_sequence_rc(&self) -> &[u8] {
        self.rev_seq_rc
    }

    fn reverse_quality_score_bytes_rc(&self) -> &[u8] {
        self.rev_quality_score_bytes_rc
    }
}

fn ensure_sequence_quality_lengths(
    mate: &'static str,
    seq_len: usize,
    qual_len: usize,
) -> Result<()> {
    if seq_len != qual_len {
        return Err(OverlapError::OrientedPairSequenceQualityLengthMismatch {
            mate,
            seq_len,
            qual_len,
        }
        .into());
    }

    Ok(())
}

#[cfg(test)]
pub(super) fn reverse_complement_bytes(seq: &[u8]) -> Box<[u8]> {
    seq.iter().rev().copied().map(complement_base).collect()
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
