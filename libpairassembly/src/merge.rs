use crate::{Read, ReadMates, ValidatedOverlap};

// Should merging be a trait?

impl<'read> ValidatedOverlap<'read> {
    /// A core, heavy-lifter function in this crate. `merge()` takes views into the original reads
    /// as well as their overlaps and generates a consensus read. It notably does not, on its own,
    /// handle quality score corrections yet. Instead, its main run is to index into the forward
    /// and reverse read mates, concatenate references to their bytes into contiguous slices, and
    /// decide which base-call should be selected for each base. Given that, this method could in
    /// principle be used for sequence read formats without quality scores like FASTA.
    pub fn merge(self) -> color_eyre::Result<Read<'read>> {
        todo!()
    }

    pub fn merge_with_correction(self) -> color_eyre::Result<Read<'read>> {
        todo!()
    }

    pub fn call_consensus(&self) -> MergeConsensus {
        let consensus = self.mates.call_consensus();
        consensus
    }
}

impl ReadMates<'_> {
    pub fn call_consensus(&self) -> MergeConsensus {
        todo!()
    }
}

#[derive(Debug)]
pub struct MergeConsensus<'overlap> {
    id: String,
    seq: &'overlap [u8],
    qual: &'overlap [u8],
}

// Should BaseOverlap instances be generated somewhere in this module?
