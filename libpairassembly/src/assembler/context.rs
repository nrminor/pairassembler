//! Internal context and overlap snapshot carriers for assembler transitions.

use std::marker::PhantomData;

use crate::{
    PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectedOrientedPair},
    merge::MergedConsensus,
    read::{OwnedReadPair, OwnedSequenceRead, ReadPair},
    validate::ValidationMetrics,
};

use super::{
    AssemblerConfig, PairInput,
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
    O: OverlapStateStorage<'pair, 'asm>,
{
    pub(super) config: &'asm AssemblerConfig,
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
pub type OverlapSearch<'asm, 'pair, R> =
    OverlapOutcome<OverlapContext<'asm, 'pair, R>, NoOverlapContext<'asm, 'pair, R>>;

/// Runtime branch produced by overlap discovery.
///
/// The fluent merge/validate/correct methods are available on the `Found` branch. Use
/// [`OverlapOutcome::and_then_found`] when no-overlap should flow through as `Ok(None)`.
#[derive(Debug, Clone)]
pub enum OverlapOutcome<Found, NoOverlap> {
    /// Overlap discovery found candidate overlap slices.
    Found(Found),
    /// Overlap discovery completed successfully but found no acceptable overlap.
    NoOverlap(NoOverlap),
}

impl<Found, NoOverlap> OverlapOutcome<Found, NoOverlap> {
    #[must_use]
    pub fn is_found(&self) -> bool {
        matches!(self, Self::Found(_))
    }

    #[must_use]
    pub fn is_no_overlap(&self) -> bool {
        matches!(self, Self::NoOverlap(_))
    }

    #[must_use]
    pub fn found(self) -> Option<Found> {
        match self {
            Self::Found(ctx) => Some(ctx),
            Self::NoOverlap(_) => None,
        }
    }

    #[must_use]
    pub fn no_overlap(self) -> Option<NoOverlap> {
        match self {
            Self::Found(_) => None,
            Self::NoOverlap(ctx) => Some(ctx),
        }
    }

    #[must_use]
    pub fn inspect_found(self, f: impl FnOnce(&Found)) -> Self {
        if let Self::Found(ctx) = &self {
            f(ctx);
        }
        self
    }

    #[must_use]
    pub fn inspect_no_overlap(self, f: impl FnOnce(&NoOverlap)) -> Self {
        if let Self::NoOverlap(ctx) = &self {
            f(ctx);
        }
        self
    }

    /// Run a fallible continuation only when overlap search found overlap slices.
    ///
    /// ```rust
    /// use libpairassembly::assembler::OverlapOutcome;
    ///
    /// # fn main() -> libpairassembly::Result<()> {
    /// let result: Option<&'static str> = OverlapOutcome::<u8, ()>::NoOverlap(())
    ///     .and_then_found(|_| Ok("not called"))?;
    ///
    /// assert_eq!(result, None);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns the error produced by the continuation. A no-overlap search result returns
    /// `Ok(None)` without invoking the continuation.
    pub fn and_then_found<T>(self, f: impl FnOnce(Found) -> Result<T>) -> Result<Option<T>> {
        match self {
            Self::Found(ctx) => f(ctx).map(Some),
            Self::NoOverlap(_) => Ok(None),
        }
    }
}

/// Merged state produced from unvalidated pair slices.
pub type MergedContext<'asm, 'pair> = MergeContext<'asm, 'pair, Unvalidated, Uncorrected>;

/// Merged state produced from validated pair slices.
pub type ValidatedMergedContext<'asm, 'pair> = MergeContext<'asm, 'pair, Validated, Uncorrected>;

/// Corrected unmerged state whose current pair slices have not been validated.
pub type CorrectedContext<'asm, 'pair, R> = CorrectedPairContext<'asm, 'pair, R, Unvalidated>;

/// Corrected unmerged state whose current pair slices have been validated.
pub type ValidatedCorrectedContext<'asm, 'pair, R> =
    CorrectedPairContext<'asm, 'pair, R, Validated>;

/// Corrected merged state produced from unvalidated slices.
pub type CorrectedMergedContext<'asm, 'pair> = CorrectedMergeContext<'asm, 'pair, Unvalidated>;

/// Corrected merged state produced from validated slices.
pub type ValidatedCorrectedMergedContext<'asm, 'pair> =
    CorrectedMergeContext<'asm, 'pair, Validated>;

/// Internal typestate carrier for merged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct MergeContext<'asm, 'pair, V, C> {
    pub(super) config: &'asm AssemblerConfig,
    pub(super) consensus: MergedConsensus,
    pub(super) overlap: PairOverlap<'pair, 'asm>,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(V, C)>,
}

/// Internal typestate carrier for corrected unmerged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct CorrectedPairContext<'asm, 'pair, R, V> {
    pub(super) config: &'asm AssemblerConfig,
    pub(super) input: &'pair PairInput<R>,
    pub(super) corrected_pair: CorrectedOrientedPair,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<V>,
}

/// Internal typestate carrier for corrected merged-stage DAG transitions.
#[derive(Debug, Clone)]
pub struct CorrectedMergeContext<'asm, 'pair, V> {
    pub(super) config: &'asm AssemblerConfig,
    pub(super) corrected_merged: CorrectedMergedRead,
    pub(super) corrected_pair: CorrectedOrientedPair,
    pub(super) validation_metrics: Option<ValidationMetrics>,
    pub(super) _marker: PhantomData<(&'pair (), V)>,
}

impl<'asm, 'pair, R, O, V, M, C> PairContext<'asm, 'pair, R, O, V, M, C>
where
    O: OverlapStateStorage<'pair, 'asm>,
{
    #[inline]
    pub(super) fn validation_metrics_ref(&self) -> Option<&ValidationMetrics> {
        self.validation_metrics.as_ref()
    }
}

impl<'asm, 'pair, R, O, V, M> PairContext<'asm, 'pair, R, O, V, M, Uncorrected>
where
    O: OverlapStateStorage<'pair, 'asm>,
{
    #[must_use]
    pub fn as_read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
    }
}

impl<'asm, 'pair, R, V, M, C> PairContext<'asm, 'pair, R, OverlapFound, V, M, C> {
    #[inline]
    #[must_use]
    pub fn overlap(&self) -> &PairOverlap<'pair, 'asm> {
        &self.overlap
    }
}

impl<'pair, R> NoOverlapContext<'_, 'pair, R> {
    #[must_use]
    pub fn read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
    }
}

impl<V, C> MergeContext<'_, '_, V, C> {
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

impl<V> CorrectedMergeContext<'_, '_, V> {
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
