//! Internal context and overlap snapshot carriers for assembler transitions.

use std::marker::PhantomData;

use crate::{PairOverlap, ReadPair, Result};

use super::{
    Assembler, PairInput,
    typestate::{HasOverlap, NoOverlap, Uncorrected, Unmerged, Unvalidated, Validated},
};

/// Internal typestate carrier for per-pair Assembler DAG transitions.
#[derive(Debug, Clone)]
pub struct PairContext<'asm, 'pair, R, O, V, M, C> {
    pub(super) assembler: &'asm Assembler,
    pub(super) input: &'pair PairInput<R>,
    pub(super) read_pair: ReadPair<'pair>,
    pub(super) overlap_outcome: OverlapOutcome,
    pub(super) _marker: PhantomData<(O, V, M, C)>,
}

#[derive(Debug, Clone, Copy)]
pub(super) enum OverlapOutcome {
    Unknown,
    Missing,
    Found(OverlapSnapshot),
}

#[derive(Debug)]
pub(super) enum OverlapBranch<C, T> {
    Value(T),
    Context(C),
}

impl<C, T> OverlapBranch<C, T> {
    pub(super) fn on_missing(self, f: impl FnOnce(C) -> Result<T>) -> Result<T> {
        match self {
            Self::Value(value) => Ok(value),
            Self::Context(ctx) => f(ctx),
        }
    }
}

/// Initial per-pair state before overlap discovery.
pub type PairReady<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, NoOverlap, Unvalidated, Unmerged, Uncorrected>;

/// Per-pair state after overlap discovery and before validation.
pub type OverlapContext<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, HasOverlap, Unvalidated, Unmerged, Uncorrected>;

/// Per-pair state after explicit overlap validation.
pub type ValidatedContext<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, HasOverlap, Validated, Unmerged, Uncorrected>;

impl<'asm, 'pair, R, O, V, M, C> PairContext<'asm, 'pair, R, O, V, M, C> {
    #[inline]
    pub(super) fn assembler_ref(&self) -> &'asm Assembler {
        self.assembler
    }

    #[inline]
    pub(super) fn read_pair_ref(&self) -> &ReadPair<'pair> {
        &self.read_pair
    }

    #[inline]
    pub(super) fn overlap_outcome(&self) -> OverlapOutcome {
        self.overlap_outcome
    }

    #[inline]
    pub(super) fn into_parts(self) -> (&'asm Assembler, &'pair PairInput<R>, ReadPair<'pair>) {
        (self.assembler, self.input, self.read_pair)
    }

    #[inline]
    pub(super) fn on_found<T>(
        self,
        f: impl FnOnce(Self, OverlapSnapshot) -> Result<T>,
    ) -> Result<OverlapBranch<Self, T>> {
        match self.overlap_outcome {
            OverlapOutcome::Found(snapshot) => Ok(OverlapBranch::Value(f(self, snapshot)?)),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => Ok(OverlapBranch::Context(self)),
        }
    }

    #[inline]
    pub(super) fn on_missing<T>(
        self,
        f: impl FnOnce(Self) -> Result<T>,
    ) -> Result<OverlapBranch<Self, T>> {
        match self.overlap_outcome {
            OverlapOutcome::Missing => Ok(OverlapBranch::Value(f(self)?)),
            OverlapOutcome::Found(_) | OverlapOutcome::Unknown => Ok(OverlapBranch::Context(self)),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct OverlapSnapshot {
    overlap_len: usize,
    r1_start_offset: usize,
    r1_end_offset: usize,
    r2_start_offset: usize,
    r2_end_offset: usize,
}

impl OverlapSnapshot {
    pub(super) fn from_overlap(overlap: &PairOverlap<'_>) -> Self {
        Self {
            overlap_len: overlap.len(),
            r1_start_offset: overlap.forward_start_offset(),
            r1_end_offset: overlap.forward_end_offset(),
            r2_start_offset: overlap.reverse_start_offset(),
            r2_end_offset: overlap.reverse_end_offset(),
        }
    }

    pub(super) fn materialize_overlap<'a>(self, pair: &'a ReadPair<'_>) -> PairOverlap<'a> {
        let r1_seq = pair.fwd_mate.sequence().as_bytes();
        let r1_qual = pair.fwd_mate.quality_scores().as_bytes();
        let r2_seq_rc = pair.rev_mate.reverse_complement();
        let mut r2_qual_rc = pair.rev_mate.quality_scores().as_bytes().to_vec();
        r2_qual_rc.reverse();

        PairOverlap::from_components(
            self.overlap_len,
            self.r1_start_offset,
            self.r1_end_offset,
            self.r2_start_offset,
            self.r2_end_offset,
            &r1_seq[self.r1_start_offset..=self.r1_end_offset],
            &r1_qual[self.r1_start_offset..=self.r1_end_offset],
            r2_seq_rc[self.r2_start_offset..=self.r2_end_offset].to_vec(),
            r2_qual_rc[self.r2_start_offset..=self.r2_end_offset].to_vec(),
        )
    }
}
