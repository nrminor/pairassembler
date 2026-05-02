//! Internal context and overlap snapshot carriers for assembler transitions.

use std::marker::PhantomData;

use crate::{
    OwnedReadPair, OwnedSequenceRead, PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedPairEvidence},
    merge::MergedConsensus,
    overlap::OverlapBounds,
    validate::ValidationMetrics,
};

use super::{
    Assembler, PairInput,
    typestate::{
        NoOverlapFound, OverlapFound, OverlapStateStorage, OverlapUnsearched, Uncorrected,
        Unmerged, Unvalidated, Validated,
    },
};

/// Internal typestate carrier for per-pair Assembler DAG transitions.
#[derive(Debug, Clone)]
#[doc(hidden)]
pub struct PairContext<'asm, 'pair, R, O, V, M, C>
where
    O: OverlapStateStorage<'pair>,
{
    pub(super) assembler: &'asm Assembler,
    pub(super) input: &'pair PairInput<R>,
    pub(super) read_pair: ReadPair<'pair>,
    pub(super) overlap: O::Storage,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(O, V, M, C)>,
}

/// Initial per-pair state before overlap discovery.
pub type PairReady<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, OverlapUnsearched, Unvalidated, Unmerged, Uncorrected>;

/// Per-pair state after overlap discovery and before validation.
pub type OverlapContext<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, OverlapFound, Unvalidated, Unmerged, Uncorrected>;

/// Per-pair state after overlap discovery found no acceptable overlap.
pub type NoOverlapContext<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, NoOverlapFound, Unvalidated, Unmerged, Uncorrected>;

/// Per-pair state after explicit overlap validation.
pub type ValidatedContext<'asm, 'pair, R> =
    PairContext<'asm, 'pair, R, OverlapFound, Validated, Unmerged, Uncorrected>;

/// Result of a successful overlap search.
#[derive(Debug, Clone)]
pub enum OverlapSearch<'asm, 'pair, R> {
    Found(OverlapContext<'asm, 'pair, R>),
    NoOverlap(NoOverlapContext<'asm, 'pair, R>),
}

impl<'asm, 'pair, R> OverlapSearch<'asm, 'pair, R> {
    #[must_use]
    pub fn is_found(&self) -> bool {
        matches!(self, Self::Found(_))
    }

    #[must_use]
    pub fn is_no_overlap(&self) -> bool {
        matches!(self, Self::NoOverlap(_))
    }

    #[must_use]
    pub fn found(self) -> Option<OverlapContext<'asm, 'pair, R>> {
        match self {
            Self::Found(ctx) => Some(ctx),
            Self::NoOverlap(_) => None,
        }
    }

    #[must_use]
    pub fn no_overlap(self) -> Option<NoOverlapContext<'asm, 'pair, R>> {
        match self {
            Self::Found(_) => None,
            Self::NoOverlap(ctx) => Some(ctx),
        }
    }

    #[must_use]
    pub fn inspect_found(self, f: impl FnOnce(&OverlapContext<'asm, 'pair, R>)) -> Self {
        if let Self::Found(ctx) = &self {
            f(ctx);
        }
        self
    }

    #[must_use]
    pub fn inspect_no_overlap(self, f: impl FnOnce(&NoOverlapContext<'asm, 'pair, R>)) -> Self {
        if let Self::NoOverlap(ctx) = &self {
            f(ctx);
        }
        self
    }

    /// Run a fallible continuation only when overlap search found overlap evidence.
    ///
    /// # Errors
    ///
    /// Returns the error produced by the continuation. A no-overlap search result returns
    /// `Ok(None)` without invoking the continuation.
    pub fn and_then_found<T>(
        self,
        f: impl FnOnce(OverlapContext<'asm, 'pair, R>) -> Result<T>,
    ) -> Result<Option<T>> {
        match self {
            Self::Found(ctx) => f(ctx).map(Some),
            Self::NoOverlap(_) => Ok(None),
        }
    }
}

/// Merged state produced from unvalidated pair evidence.
pub type MergedContext<'asm, 'pair> = MergeContext<'asm, 'pair, Unvalidated, Uncorrected>;

/// Merged state produced from validated pair evidence.
pub type ValidatedMergedContext<'asm, 'pair> = MergeContext<'asm, 'pair, Validated, Uncorrected>;

/// Corrected unmerged state whose current pair evidence has not been validated.
pub type CorrectedContext<'asm, 'pair, R> = CorrectedPairContext<'asm, 'pair, R, Unvalidated>;

/// Corrected unmerged state whose current pair evidence has been validated.
pub type ValidatedCorrectedContext<'asm, 'pair, R> =
    CorrectedPairContext<'asm, 'pair, R, Validated>;

/// Corrected merged state produced from unvalidated evidence.
pub type CorrectedMergedContext<'asm> = CorrectedMergeContext<'asm, Unvalidated>;

/// Corrected merged state produced from validated evidence.
pub type ValidatedCorrectedMergedContext<'asm> = CorrectedMergeContext<'asm, Validated>;

/// Internal typestate carrier for merged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct MergeContext<'asm, 'pair, V, C> {
    pub(super) assembler: &'asm Assembler,
    pub(super) consensus: MergedConsensus,
    pub(super) overlap: PairOverlap<'pair>,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(V, C)>,
}

/// Internal typestate carrier for corrected unmerged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct CorrectedPairContext<'asm, 'pair, R, V> {
    pub(super) assembler: &'asm Assembler,
    pub(super) input: &'pair PairInput<R>,
    pub(super) corrected_pair: CorrectedPairEvidence,
    pub(super) overlap_bounds: OverlapBounds,
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

impl<'pair, R, O, V, M, C> PairContext<'_, 'pair, R, O, V, M, C>
where
    O: OverlapStateStorage<'pair>,
{
    #[inline]
    pub(super) fn read_pair_ref(&self) -> &ReadPair<'pair> {
        &self.read_pair
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }
}

impl<'pair, R, O, V, M> PairContext<'_, 'pair, R, O, V, M, Uncorrected>
where
    O: OverlapStateStorage<'pair>,
{
    #[must_use]
    pub fn as_read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
    }
}

impl<'pair, R, V, M, C> PairContext<'_, 'pair, R, OverlapFound, V, M, C> {
    #[inline]
    #[must_use]
    pub fn overlap(&self) -> &PairOverlap<'pair> {
        &self.overlap
    }
}

impl<'pair, R> NoOverlapContext<'_, 'pair, R> {
    #[must_use]
    pub fn read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
    }
}

impl<'pair, V, C> MergeContext<'_, 'pair, V, C> {
    #[inline]
    pub(super) fn consensus_ref(&self) -> &MergedConsensus {
        &self.consensus
    }

    #[inline]
    pub(super) fn overlap_ref(&self) -> &PairOverlap<'pair> {
        &self.overlap
    }

    #[must_use]
    pub fn id(&self) -> &str {
        self.consensus.id()
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        self.consensus.sequence()
    }

    #[must_use]
    pub fn quality_score_bytes(&self) -> &[u8] {
        self.consensus.quality_score_bytes()
    }

    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }
}

impl<V> MergeContext<'_, '_, V, Uncorrected> {
    /// Consume this merged context into an owned FASTQ-shaped read.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_read(self) -> Result<OwnedSequenceRead> {
        OwnedSequenceRead::try_from(self.consensus)
    }
}

impl<R, V> CorrectedPairContext<'_, '_, R, V> {
    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }

    /// Consume this corrected pair context into an owned FASTQ-shaped read pair.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_pair(self) -> Result<OwnedReadPair> {
        OwnedReadPair::try_from(self.corrected_pair)
    }
}

impl<V> CorrectedMergeContext<'_, V> {
    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }

    /// Consume this corrected merged context into an owned FASTQ-shaped read.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_read(self) -> Result<OwnedSequenceRead> {
        OwnedSequenceRead::try_from(self.corrected_merged)
    }
}
