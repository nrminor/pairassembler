/// Public convenience re-exports for common `libpairassembly` use cases.
pub use crate::{
    Error,
    assembler::{
        Assembler, AssemblerBuilder, AssemblerConfig, MergeParams, MergeTiePolicy, PairInput,
        PairReady, SeqRecordView,
    },
    correct::{CorrectedMergedRead, CorrectionParams},
    errors::Result,
    merge::MergedRead,
    overlap::{OverlapParams, PairOverlap, TiePolicy},
    read::{OwnedReadPair, OwnedSequenceRead, ReadPair, SequenceRead},
    validate::{OverlapValidator, ValidatedOverlap},
};

#[macro_use]
pub mod utils {
    pub(crate) const PHRED_OFFSET: u8 = 33;

    #[inline]
    pub(crate) fn fastq_ascii_to_phred(quality: u8) -> u8 {
        quality.saturating_sub(PHRED_OFFSET)
    }

    #[inline]
    pub(crate) fn phred_to_fastq_ascii(phred: u8) -> u8 {
        phred.saturating_add(PHRED_OFFSET)
    }

    pub(crate) fn decode_fastq_quality_scores(qualities: &[u8]) -> Box<[u8]> {
        qualities
            .iter()
            .copied()
            .map(fastq_ascii_to_phred)
            .collect()
    }

    pub(crate) fn encode_fastq_quality_scores_in_place(qualities: &mut [u8]) {
        for quality in qualities {
            *quality = phred_to_fastq_ascii(*quality);
        }
    }

    /// Compute the reverse complement of a DNA sequence, preserving case and supporting all IUPAC bases.
    ///
    /// # Panics
    ///
    /// Panics when `seq` contains a non-IUPAC DNA base.
    pub fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                'a' => 't',
                't' => 'a',
                'c' => 'g',
                'g' => 'c',
                'R' => 'Y',
                'Y' => 'R',
                'S' => 'S',
                'W' => 'W',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'D' => 'H',
                'H' => 'D',
                'V' => 'B',
                'N' => 'N',
                'r' => 'y',
                'y' => 'r',
                's' => 's',
                'w' => 'w',
                'k' => 'm',
                'm' => 'k',
                'b' => 'v',
                'd' => 'h',
                'h' => 'd',
                'v' => 'b',
                'n' => 'n',
                invalid => panic!("Invalid DNA base encountered in sequence: '{invalid}'"),
            })
            .collect()
    }

    /// Construct a [`SequenceRead`](crate::read::SequenceRead) from string literals after checking
    /// at compile time that sequence and quality literals have matching lengths.
    #[macro_export]
    macro_rules! new_sequence_read {
        ($id:expr, $seq:expr, $qual:expr) => {{
            const _: () = {
                ["quality and sequence length mismatch"]
                    [(($seq.len() == $qual.len()) as usize) ^ 1];
            };
            $crate::read::SequenceRead::from_literal_parts($id, $seq, $qual)
        }};
    }
}
pub use utils::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_handles_iupac() {
        let seq = "ACGTRYSWKMBDHVN";
        let rc = reverse_complement(seq);
        assert_eq!(rc, "NBDHVKMWSRYACGT");
    }
}
