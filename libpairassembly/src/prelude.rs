/// `prelude` is first and foremost for re-exporting the core parts of the rest of the library's
/// modules, and also for core types shared across those libraries.
//
// RE-EXPORTS
// ------------------------------------------------------------------------------------------------
pub use crate::{
    Error,
    errors::Result,
    overlap::{MateOverlap, OverlapParams},
    validate::{BaseCallValidator, ValidatedOverlap},
};
// ------------------------------------------------------------------------------------------------

// OPTIONAL RE-EXPORTS
// ------------------------------------------------------------------------------------------------
#[cfg(feature = "noodles")]
pub use crate::io::noodles::*;
// #[cfg(feature = "rust-bio")]
// pub use crate::io::rust_bio::*;
// #[cfg(feature = "needletail")]
// pub use crate::io::needletail::*;
// #[cfg(feature = "binseq")]
// pub use crate::io::binseq::*;
// ------------------------------------------------------------------------------------------------

// CORE TYPES
// ------------------------------------------------------------------------------------------------

use crate::errors::{PairingError, SequenceQualityLengthMismatch};

#[derive(Debug)]
pub struct SequenceRead<'read> {
    id: &'read str,
    seq: &'read str,
    qual: &'read str,
}

impl<'read> SequenceRead<'read> {
    pub(crate) fn new(id: &'read str, seq: &'read str, qual: &'read str) -> Self {
        assert_eq!(seq.len(), qual.len());
        SequenceRead { id, seq, qual }
    }

    pub fn try_new(id: &'read str, seq: &'read str, qual: &'read str) -> Result<Self> {
        if seq.len() != qual.len() {
            return Err(SequenceQualityLengthMismatch(
                seq.to_string(),
                seq.len(),
                qual.to_string(),
                qual.len(),
            )
            .into());
        }

        Ok(SequenceRead { id, seq, qual })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        assert_eq!(self.seq.len(), self.qual.len());
        self.seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn reverse_complement(&self) -> Vec<u8> {
        let rc = reverse_complement(self.seq);
        rc.as_bytes().to_vec()
    }

    #[inline]
    #[must_use]
    pub fn check_for_mate(&self, possible_mate: &SequenceRead) -> bool {
        self.id == possible_mate.id
    }

    #[inline]
    #[must_use]
    pub fn id(&self) -> &str {
        self.id
    }

    #[inline]
    #[must_use]
    pub fn quality_scores(&self) -> &str {
        self.qual
    }

    #[inline]
    #[must_use]
    pub fn sequence(&self) -> &str {
        self.seq
    }
}

#[derive(Debug)]
pub struct ReadMates<'mate> {
    pub fwd_mate: SequenceRead<'mate>,
    pub rev_mate: SequenceRead<'mate>,
}

impl<'a> ReadMates<'a> {
    pub fn from(read1: SequenceRead<'a>, read2: SequenceRead<'a>) -> Result<Self> {
        if read1.id != read2.id {
            return Err(
                PairingError::UnmatchedIds(read1.id.to_string(), read2.id.to_string()).into(),
            );
        }
        let pair = ReadMates {
            fwd_mate: read1,
            rev_mate: read2,
        };
        Ok(pair)
    }
}

#[macro_use]
pub mod utils {

    // TODO: refactor so that String and Vec heap allocations don't need to be performed as redundantly.
    /// Compute the reverse complement of a DNA sequence, preserving case and supporting all IUPAC bases.
    /// Panics on invalid input.
    ///
    pub fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                // Uppercase unambiguous
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',

                // Lowercase unambiguous
                'a' => 't',
                't' => 'a',
                'c' => 'g',
                'g' => 'c',

                // Uppercase ambiguous
                'R' => 'Y', // A or G -> T or C
                'Y' => 'R', // C or T -> G or A
                'S' => 'S', // G or C -> C or G (self-complementary)
                'W' => 'W', // A or T -> T or A (self-complementary)
                'K' => 'M', // G or T -> C or A
                'M' => 'K', // A or C -> T or G
                'B' => 'V', // C or G or T -> G or C or A
                'D' => 'H', // A or G or T -> T or C or A
                'H' => 'D', // A or C or T -> T or G or A
                'V' => 'B', // A or C or G -> T or G or C
                'N' => 'N', // any -> any

                // Lowercase ambiguous
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

    /// Macro that should be used to construct sequence read instances with string literals known
    /// at compile-time. The macro mainly exists to ensure that the sequence literals for the sequence
    /// and the quality scores are the same length; if they are not, it will prevent the user from
    /// constructing the instance.
    #[macro_export]
    macro_rules! new_sequence_read {
        ($id:expr, $seq:expr, $qual:expr) => {{
            const _: () = {
                ["quality and sequence length mismatch"]
                    [(($seq.len() == $qual.len()) as usize) ^ 1];
            };
            SequenceRead {
                id: $id,
                seq: $seq,
                qual: $qual,
            }
        }};
    }
}
pub use utils::*;

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    #[test]
    fn test_reverse_complement_handles_iupac() {
        let seq = "ACGTRYSWKMBDHVN";
        let rc = reverse_complement(seq);
        assert_eq!(rc, "NBDHVKMWSRYACGT");
    }

    #[test]
    fn test_try_new_rejects_seq_qual_length_mismatch() {
        let result = SequenceRead::try_new("r1", "ACGT", "III");
        assert!(result.is_err());
    }

    #[test]
    fn test_read_mates_from_rejects_mismatched_ids() {
        let r1 = SequenceRead::new("read1", "ACGT", "IIII");
        let r2 = SequenceRead::new("read2", "ACGT", "IIII");
        let result = ReadMates::from(r1, r2);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_mates_from_accepts_matching_ids() {
        let r1 = SequenceRead::new("read1", "ACGT", "IIII");
        let r2 = SequenceRead::new("read1", "TGCA", "IIII");
        let result = ReadMates::from(r1, r2);
        assert!(result.is_ok());
    }
}
