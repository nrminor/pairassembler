use crate::{
    Result, errors::OverlapError, prelude::utils::decode_fastq_quality_scores, read::ReadPair,
};

use super::{OverlapBounds, OverlapFinder, PairOverlap, private};

impl<'a> ReadPair<'a> {
    /// Discover the best overlap between this read pair.
    ///
    /// Returns `Ok(None)` when no overlap candidate satisfies the configured thresholds.
    ///
    /// # Errors
    ///
    /// Returns an error when overlap candidate reconciliation fails (for example, tie rejection),
    /// or when computed overlap coordinates are inconsistent with read bounds.
    pub fn overlap(&self, params: &super::OverlapParams) -> Result<Option<PairOverlap<'a>>> {
        OverlapFinder::new(params).find(*self)
    }

    pub(crate) fn to_oriented_slices(self) -> OrientedPairSlices<'a> {
        OrientedPairSlices::from_fastq_ascii_parts(
            self.fwd_id(),
            self.fwd_sequence_bytes(),
            self.fwd_quality_bytes(),
            self.rev_sequence_bytes(),
            self.rev_quality_bytes(),
        )
    }
}

impl<'a> OrientedPairSlices<'a> {
    fn from_fastq_ascii_parts(
        id: &'a str,
        fwd_seq: &'a [u8],
        fwd_qual: &[u8],
        rev_raw_seq: &'a [u8],
        rev_raw_qual: &[u8],
    ) -> Self {
        let fwd_qual = decode_fastq_quality_scores(fwd_qual);
        let mut rev_qual_rev = decode_fastq_quality_scores(rev_raw_qual);
        rev_qual_rev.reverse();
        let rev_seq_rc = reverse_complement_bytes(rev_raw_seq);

        Self {
            id,
            fwd_seq,
            fwd_qual,
            rev_seq_rc,
            rev_qual_rev,
        }
    }
}

#[derive(Debug, Clone)]
/// Retained pair slices after converting reads into overlap orientation.
///
/// Forward sequence bytes are borrowed from the original read. Quality scores and reverse-mate
/// orientation are owned fixed-size buffers because FASTQ ASCII qualities must be decoded and the
/// reverse mate must be reverse-complemented before downstream validation, merge, and correction.
pub(crate) struct OrientedPairSlices<'a> {
    pub(crate) id: &'a str,
    pub(crate) fwd_seq: &'a [u8],
    pub(crate) fwd_qual: Box<[u8]>,
    pub(crate) rev_seq_rc: Box<[u8]>,
    pub(crate) rev_qual_rev: Box<[u8]>,
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

impl private::Sealed for OrientedPairSlices<'_> {}

impl HasOrientedPairSlices for OrientedPairSlices<'_> {
    fn pair_id(&self) -> &str {
        self.id
    }

    fn forward_sequence(&self) -> &[u8] {
        self.fwd_seq
    }

    fn forward_quality_score_bytes(&self) -> &[u8] {
        &self.fwd_qual
    }

    fn reverse_sequence_rc(&self) -> &[u8] {
        &self.rev_seq_rc
    }

    fn reverse_quality_score_bytes_rc(&self) -> &[u8] {
        &self.rev_qual_rev
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

pub(super) fn reverse_complement_bytes(seq: &[u8]) -> Box<[u8]> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}
