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

impl<'pair, R, O, V, M, C> PairContext<'_, 'pair, R, O, V, M, C> {
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
}

impl<'pair, R, O, V, M> PairContext<'_, 'pair, R, O, V, M, Uncorrected> {
    #[must_use]
    pub fn as_read_pair(&self) -> ReadPair<'pair> {
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
