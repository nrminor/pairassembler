//! Iterator adaptor for processing paired records through an [`Assembler`].

use crate::{Result, read::OwnedSequenceRead};

use super::{Assembler, PairInput, SeqRecordView};

#[derive(Debug)]
pub struct ProcessIter<'asm, I> {
    pub(super) assembler: &'asm Assembler,
    pub(super) iter: I,
}

impl<I, R> Iterator for ProcessIter<'_, I>
where
    I: Iterator<Item = PairInput<R>>,
    R: SeqRecordView,
{
    type Item = Result<Option<OwnedSequenceRead>>;

    fn next(&mut self) -> Option<Self::Item> {
        let pair = self.iter.next()?;
        Some(self.assembler.process_pair(&pair))
    }
}
