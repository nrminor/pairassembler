pub use crate::validate::ValidatedOverlap;

#[derive(Debug)]
pub struct Read<'read> {
    id: &'read str,
    seq: &'read str,
    qual: &'read str,
}

impl<'read> Read<'read> {
    pub fn new(id: &'read str, seq: &'read str, qual: &'read str) -> Self {
        assert_eq!(seq.len(), qual.len());
        Self { id, seq, qual }
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.seq.len(), self.qual.len());
        self.seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn reverse_complement(&self) -> String {
        todo!()
    }

    #[inline]
    pub fn check_for_mate(&self, possible_mate: &Read) -> bool {
        self.id == possible_mate.id
    }

    #[inline]
    pub fn id(&self) -> &str {
        self.id
    }

    #[inline]
    pub fn quality_scores(&self) -> &str {
        self.qual
    }

    #[inline]
    pub fn sequence(&self) -> &str {
        self.seq
    }
}

#[derive(Debug)]
pub struct ReadMates<'mate> {
    pub fwd_mate: Read<'mate>,
    pub rev_mate: Read<'mate>,
}

impl ReadMates<'_> {
    pub fn new(read1: &Read, read2: &Read) -> color_eyre::Result<Self> {
        todo!()
    }
}

pub mod utils {
    /// Compute the reverse complement of a DNA sequence.
    pub fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' | 'a' => 'T',
                'T' | 't' => 'A',
                'C' | 'c' => 'G',
                'G' | 'g' => 'C',
                other => todo!(),
            })
            .collect()
    }
}
