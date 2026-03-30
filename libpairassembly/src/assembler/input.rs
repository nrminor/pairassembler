//! Pair input adapters for assembler entrypoints.

use crate::{ReadPair, Result, SequenceRead};

use super::SeqRecordView;

/// Pair wrapper accepted by assembler entrypoints.
#[derive(Debug)]
pub struct PairInput<R> {
    pub r1: R,
    pub r2: R,
}

impl<R> PairInput<R> {
    /// Construct a paired input wrapper.
    #[must_use]
    pub fn new(r1: R, r2: R) -> Self {
        Self { r1, r2 }
    }

    /// Convert a generic pair input into the canonical internal [`ReadPair`] form.
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
