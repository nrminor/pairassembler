//! Internal context and overlap snapshot carriers for assembler transitions.

use std::marker::PhantomData;

use crate::{
    OwnedReadPair, OwnedSequenceRead, PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedReadPair},
    merge::MergedRead,
    validate::ValidationMetrics,
};

use super::{
    Assembler, IntoOwnedPairRecordParts, IntoOwnedRecordParts, PairInput,
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
    Found(PairOverlap<'pair>),
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

/// Merged state produced from unvalidated pair evidence.
pub type MergedContext<'asm> = MergeContext<'asm, Unvalidated, Uncorrected>;

/// Merged state produced from validated pair evidence.
pub type ValidatedMergedContext<'asm> = MergeContext<'asm, Validated, Uncorrected>;

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
        f: impl FnOnce(Self, PairOverlap<'pair>) -> Result<T>,
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

impl<'asm, 'pair, R, O, V, M> PairContext<'asm, 'pair, R, O, V, M, Uncorrected> {
    #[must_use]
    pub fn as_read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
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
    pub fn quality_score_bytes(&self) -> &[u8] {
        self.merged.quality_score_bytes()
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

impl<'asm, V> MergeContext<'asm, V, Uncorrected> {
    /// Consume this merged context into an owned FASTQ-shaped read.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_read(self) -> Result<OwnedSequenceRead> {
        let (id, seq, qual) = self.merged.into_owned_record_parts();
        OwnedSequenceRead::try_from_parts(id, seq, qual)
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

    /// Consume this corrected pair context into an owned FASTQ-shaped read pair.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_pair(self) -> Result<OwnedReadPair> {
        let (id, fwd_seq, fwd_qual, rev_seq, rev_qual) =
            self.corrected_pair.into_owned_pair_record_parts();
        OwnedReadPair::try_from_parts(id, fwd_seq, fwd_qual, rev_seq, rev_qual)
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

    /// Consume this corrected merged context into an owned FASTQ-shaped read.
    ///
    /// # Errors
    ///
    /// Returns an error if an internal invariant is violated and sequence or quality bytes are not
    /// valid UTF-8.
    pub fn into_owned_read(self) -> Result<OwnedSequenceRead> {
        let (id, seq, qual) = self.corrected_merged.into_owned_record_parts();
        OwnedSequenceRead::try_from_parts(id, seq, qual)
    }
}
