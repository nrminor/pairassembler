//! Pair input adapters for assembler entrypoints.

use crate::{
    Result,
    read::{ReadPair, SequenceRead},
};

use super::SeqRecordView;

/// Pair wrapper accepted by assembler entrypoints.
///
/// `PairInput` keeps parser integration at the application boundary: each mate can be any type that
/// implements [`SeqRecordView`]. Conversion into the library's canonical borrowed [`ReadPair`] is
/// checked when processing begins.
#[derive(Debug)]
pub struct PairInput<R> {
    /// First mate, conventionally R1/forward.
    pub r1: R,
    /// Second mate, conventionally R2/reverse.
    pub r2: R,
}

impl<R> PairInput<R> {
    /// Construct a paired input wrapper.
    ///
    /// ```rust
    /// use libpairassembly::prelude::*;
    ///
    /// # fn main() -> libpairassembly::Result<()> {
    /// let pair = PairInput::new(
    ///     SequenceRead::try_new("read-1", "ACGT", "IIII")?,
    ///     SequenceRead::try_new("read-1", "TGCA", "IIII")?,
    /// );
    /// assert_eq!(pair.try_into_read_pair()?.fwd_id(), "read-1");
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn new(r1: R, r2: R) -> Self {
        Self { r1, r2 }
    }

    /// Convert a generic pair input into the canonical internal [`ReadPair`] form.
    ///
    /// This checks sequence/quality length equality for each mate and then checks that the two mate
    /// identifiers match.
    ///
    /// # Errors
    ///
    /// Returns an error if either record has invalid sequence/quality structure or
    /// if IDs do not correspond to a valid read pair.
    pub fn try_into_read_pair(&self) -> Result<ReadPair<'_>>
    where
        R: SeqRecordView,
    {
        let read1 = SequenceRead::try_new(self.r1.id(), self.r1.seq(), self.r1.qual())?;
        let read2 = SequenceRead::try_new(self.r2.id(), self.r2.seq(), self.r2.qual())?;
        ReadPair::from(read1, read2)
    }
}
