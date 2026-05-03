//! Boundary traits for assembler input interop.

use crate::read::SequenceRead;

/// Boundary trait for pair records accepted by the assembler API.
///
/// Implement this for parser-owned or application-owned record types to use [`Assembler`] without
/// first copying them into [`SequenceRead`]. The quality string returned by [`SeqRecordView::qual`]
/// must be FASTQ ASCII text with the same length as the sequence string.
///
/// [`Assembler`]: crate::assembler::Assembler
///
/// ```rust
/// use libpairassembly::assembler::{PairInput, SeqRecordView};
///
/// struct ParserRecord<'a> {
///     id: &'a str,
///     seq: &'a str,
///     qual: &'a str,
/// }
///
/// impl SeqRecordView for ParserRecord<'_> {
///     fn id(&self) -> &str {
///         self.id
///     }
///
///     fn seq(&self) -> &str {
///         self.seq
///     }
///
///     fn qual(&self) -> &str {
///         self.qual
///     }
/// }
///
/// # fn main() -> libpairassembly::Result<()> {
/// let pair = PairInput::new(
///     ParserRecord { id: "read-1", seq: "ACGT", qual: "IIII" },
///     ParserRecord { id: "read-1", seq: "TGCA", qual: "IIII" },
/// );
/// let read_pair = pair.try_into_read_pair()?;
/// assert_eq!(read_pair.fwd_sequence(), "ACGT");
/// # Ok(())
/// # }
/// ```
pub trait SeqRecordView {
    /// Return the pair identifier for this record.
    fn id(&self) -> &str;

    /// Return the nucleotide sequence as text.
    fn seq(&self) -> &str;

    /// Return the FASTQ ASCII quality string.
    fn qual(&self) -> &str;
}

impl SeqRecordView for SequenceRead<'_> {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &str {
        self.sequence()
    }

    fn qual(&self) -> &str {
        self.quality_scores()
    }
}
