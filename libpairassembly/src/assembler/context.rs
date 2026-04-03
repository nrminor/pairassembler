//! Internal context and overlap snapshot carriers for assembler transitions.

use std::marker::PhantomData;

use crate::{
    PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedReadPair},
    merge::MergedRead,
    overlap::PreparedPair,
    validate::ValidationMetrics,
};

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
    pub(super) overlap_outcome: OverlapOutcome<'pair>,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(O, V, M, C)>,
}

#[derive(Debug, Clone)]
pub(super) enum OverlapOutcome<'pair> {
    Unknown,
    Missing,
    Found(FoundOverlap<'pair>),
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

/// Merged state after unchecked merge.
pub type MergedContext<'asm> = MergeContext<'asm, Unvalidated, Uncorrected>;

/// Merged state after validation-aware merge.
pub type ValidatedMergedContext<'asm> = MergeContext<'asm, Validated, Uncorrected>;

/// Corrected unmerged state after correction without prior validation.
pub type CorrectedContext<'asm, 'pair, R> = CorrectedPairContext<'asm, 'pair, R, Unvalidated>;

/// Corrected unmerged state after correction with prior validation.
pub type ValidatedCorrectedContext<'asm, 'pair, R> =
    CorrectedPairContext<'asm, 'pair, R, Validated>;

/// Corrected merged state after correction without prior validation.
pub type CorrectedMergedContext<'asm> = CorrectedMergeContext<'asm, Unvalidated>;

/// Corrected merged state after correction with prior validation.
pub type ValidatedCorrectedMergedContext<'asm> = CorrectedMergeContext<'asm, Validated>;

/// Internal typestate carrier for merged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct MergeContext<'asm, V, C> {
    pub(super) assembler: &'asm Assembler,
    pub(super) merged: MergedRead,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(V, C)>,
}

/// Internal typestate carrier for corrected unmerged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct CorrectedPairContext<'asm, 'pair, R, V> {
    pub(super) assembler: &'asm Assembler,
    pub(super) input: &'pair PairInput<R>,
    pub(super) corrected_pair: CorrectedReadPair,
    pub(super) overlap_outcome: OverlapOutcome<'pair>,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<V>,
}

/// Internal typestate carrier for corrected merged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct CorrectedMergeContext<'asm, V> {
    pub(super) assembler: &'asm Assembler,
    pub(super) corrected_merged: CorrectedMergedRead,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<V>,
}

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
    pub(super) fn overlap_outcome(&self) -> &OverlapOutcome<'pair> {
        &self.overlap_outcome
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }

    #[inline]
    pub(super) fn into_parts(
        self,
    ) -> (
        &'asm Assembler,
        &'pair PairInput<R>,
        ReadPair<'pair>,
        Option<ValidationMetrics>,
    ) {
        (
            self.assembler,
            self.input,
            self.read_pair,
            self.validation_metrics,
        )
    }

    #[inline]
    pub(super) fn on_found<T>(
        self,
        f: impl FnOnce(Self, FoundOverlap<'pair>) -> Result<T>,
    ) -> Result<OverlapBranch<Self, T>> {
        match self.overlap_outcome.clone() {
            OverlapOutcome::Found(found) => Ok(OverlapBranch::Value(f(self, found)?)),
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

impl<'asm, V, C> MergeContext<'asm, V, C> {
    #[inline]
    pub(super) fn assembler_ref(&self) -> &'asm Assembler {
        self.assembler
    }

    #[inline]
    pub(super) fn merged_ref(&self) -> &MergedRead {
        &self.merged
    }

    #[must_use]
    pub fn id(&self) -> &str {
        self.merged.id()
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        self.merged.sequence()
    }

    #[must_use]
    pub fn qualities(&self) -> &[u8] {
        self.merged.qualities()
    }

    #[inline]
    pub(super) fn into_parts(self) -> (&'asm Assembler, MergedRead, Option<ValidationMetrics>) {
        (self.assembler, self.merged, self.validation_metrics)
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }
}

impl<'asm, 'pair, R, V> CorrectedPairContext<'asm, 'pair, R, V> {
    #[inline]
    pub(super) fn assembler_ref(&self) -> &'asm Assembler {
        self.assembler
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }

    pub fn into_corrected_read_pair(self) -> CorrectedReadPair {
        self.corrected_pair
    }
}

impl<'asm, V> CorrectedMergeContext<'asm, V> {
    #[inline]
    pub(super) fn assembler_ref(&self) -> &'asm Assembler {
        self.assembler
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }

    pub fn into_corrected_merged_read(self) -> CorrectedMergedRead {
        self.corrected_merged
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct OverlapBounds {
    overlap_len: usize,
    r1_start_offset: usize,
    r2_start_offset: usize,
}

#[derive(Debug, Clone)]
pub(super) struct FoundOverlap<'pair> {
    prepared: PreparedPair<'pair>,
    bounds: OverlapBounds,
}

impl OverlapBounds {
    pub(super) fn from_overlap(overlap: &PairOverlap<'_>) -> Self {
        Self {
            overlap_len: overlap.len(),
            r1_start_offset: overlap.forward_start_offset(),
            r2_start_offset: overlap.reverse_start_offset(),
        }
    }

    #[inline]
    pub(super) fn overlap_len(self) -> usize {
        self.overlap_len
    }

    #[inline]
    pub(super) fn fwd_start_offset(self) -> usize {
        self.r1_start_offset
    }

    #[inline]
    pub(super) fn fwd_end_offset(self) -> usize {
        self.r1_start_offset + self.overlap_len - 1
    }

    #[inline]
    pub(super) fn rev_start_offset(self) -> usize {
        self.r2_start_offset
    }

    #[inline]
    pub(super) fn rev_end_offset(self) -> usize {
        self.r2_start_offset + self.overlap_len - 1
    }
}

impl<'pair> FoundOverlap<'pair> {
    pub(super) fn new(prepared: PreparedPair<'pair>, bounds: OverlapBounds) -> Self {
        Self { prepared, bounds }
    }

    pub(super) fn materialize_overlap(&self) -> PairOverlap<'_> {
        let bounds = self.bounds;

        PairOverlap::try_new(
            bounds.overlap_len(),
            bounds.fwd_start_offset(),
            bounds.fwd_end_offset(),
            bounds.rev_start_offset(),
            bounds.rev_end_offset(),
            &self.prepared.fwd_seq[bounds.fwd_start_offset()..=bounds.fwd_end_offset()],
            &self.prepared.fwd_qual[bounds.fwd_start_offset()..=bounds.fwd_end_offset()],
            self.prepared.rev_seq_rc[bounds.rev_start_offset()..=bounds.rev_end_offset()].to_vec(),
            self.prepared.rev_qual_rev[bounds.rev_start_offset()..=bounds.rev_end_offset()]
                .to_vec(),
        )
        .expect("retained found overlaps should always materialize into valid overlap windows")
    }

    pub(super) fn prepared_pair(&self) -> &PreparedPair<'pair> {
        &self.prepared
    }

    pub(super) fn bounds(&self) -> OverlapBounds {
        self.bounds
    }
}
