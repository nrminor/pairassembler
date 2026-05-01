//! Iterator adaptor for processing paired records through an [`Assembler`].

use crate::{OwnedSequenceRead, Result};

use super::{Assembler, ExecutionPolicy, PairInput, SeqRecordView};

#[derive(Debug)]
pub struct ProcessIter<'asm, I> {
    pub(super) assembler: &'asm Assembler,
    pub(super) iter: I,
    pub(super) execution: ExecutionPolicy,
}

impl<I, R> Iterator for ProcessIter<'_, I>
where
    I: Iterator<Item = PairInput<R>>,
    R: SeqRecordView,
{
    type Item = Result<OwnedSequenceRead>;

    fn next(&mut self) -> Option<Self::Item> {
        let pair = self.iter.next()?;
        let result = match self.execution {
            ExecutionPolicy::Record => self.assembler.process_pair(&pair),
            ExecutionPolicy::Batch(_policy) => self.assembler.process_pair(&pair),
        };
        Some(result)
    }
}
