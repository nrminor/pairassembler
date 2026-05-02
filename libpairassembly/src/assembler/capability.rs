//! Internal capability traits for post-overlap operation contracts.

use crate::{
    PairOverlap, ReadPair, Result,
    correct::{CorrectedMergedRead, CorrectedReadPair, CorrectionWindow},
    errors::OverlapError,
    merge::{MergeView, MergedRead},
    overlap::{OverlapBounds, PreparedPair},
    validate::{ValidatedOverlap, ValidationMetrics},
};

use super::{
    CorrectedMergeContext, CorrectedPairContext, MergeContext, PairContext, context::OverlapOutcome,
};

pub(crate) mod private {
    pub(crate) trait Sealed {}
}

/// Internal marker trait for state/output carriers participating in assembler DAG operations.
pub(crate) trait PairState: private::Sealed {}

/// Capability for producing canonical overlap evidence for the current pair state.
pub(crate) trait HasPairOverlap: PairState {
    fn pair_overlap(&self) -> Result<PairOverlap<'_>>;
}

/// Capability for borrowing source read-pair evidence.
pub(crate) trait HasReadPair: PairState {
    fn read_pair(&self) -> ReadPair<'_>;
}

/// Capability for exposing normalized merge-ready overlap views.
pub(crate) trait HasMergeableOverlap: HasReadPair + HasPairOverlap {
    fn merge_view(&self) -> Result<MergeView<'_>>;
}

impl<R, O, V, M, C> HasMergeableOverlap for PairContext<'_, '_, R, O, V, M, C> {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        let overlap = match self.overlap_outcome() {
            OverlapOutcome::Found(overlap) => overlap,
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                return Err(OverlapError::NoOverlapFound.into());
            },
        };

        MergeView::from_pair_overlap(overlap)
    }
}

impl<R, V> HasMergeableOverlap for CorrectedPairContext<'_, '_, R, V> {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        let pair = self.read_pair();
        let overlap = match &self.overlap_outcome {
            OverlapOutcome::Found(overlap) => overlap,
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                return Err(OverlapError::NoOverlapFound.into());
            },
        };

        merge_view_from_fastq_pair_bounds(pair, overlap.bounds())
    }
}

impl HasMergeableOverlap for ValidatedOverlap<'_> {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        MergeView::from_pair_overlap(self.overlap())
    }
}

fn merge_view_from_fastq_pair_bounds(
    pair: ReadPair<'_>,
    bounds: OverlapBounds,
) -> Result<MergeView<'_>> {
    MergeView::from_pair_and_bounds(pair, bounds)
}

/// Capability for exposing an aligned overlap-local correction window.
pub(crate) trait HasCorrectionWindow: PairState {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>>;
}

/// Capability for exposing consensus record payload.
pub(crate) trait HasConsensusRecord: PairState {
    fn consensus_id(&self) -> &str;
    fn consensus_seq(&self) -> &[u8];
    fn consensus_quality_score_bytes(&self) -> &[u8];
}

/// Capability for exposing retained validation-stage metrics.
pub(crate) trait HasValidationMetrics: PairState {
    fn validation_metrics(&self) -> &ValidationMetrics;
}

impl<R, O, V, M, C> private::Sealed for PairContext<'_, '_, R, O, V, M, C> {}
impl<R, O, V, M, C> PairState for PairContext<'_, '_, R, O, V, M, C> {}

impl private::Sealed for ValidatedOverlap<'_> {}
impl PairState for ValidatedOverlap<'_> {}

impl private::Sealed for MergedRead {}
impl PairState for MergedRead {}

impl<V, C> private::Sealed for MergeContext<'_, V, C> {}
impl<V, C> PairState for MergeContext<'_, V, C> {}

impl<R, V> private::Sealed for CorrectedPairContext<'_, '_, R, V> {}
impl<R, V> PairState for CorrectedPairContext<'_, '_, R, V> {}

impl<V> private::Sealed for CorrectedMergeContext<'_, V> {}
impl<V> PairState for CorrectedMergeContext<'_, V> {}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl private::Sealed for CorrectedReadPair {}
impl PairState for CorrectedReadPair {}

impl<R, O, V, M, C> HasReadPair for PairContext<'_, '_, R, O, V, M, C> {
    fn read_pair(&self) -> ReadPair<'_> {
        *self.read_pair_ref()
    }
}

impl<R, O, V, M, C> HasPairOverlap for PairContext<'_, '_, R, O, V, M, C> {
    fn pair_overlap(&self) -> Result<PairOverlap<'_>> {
        match self.overlap_outcome() {
            OverlapOutcome::Found(overlap) => Ok(overlap.clone()),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }
}

impl<R, V> HasReadPair for CorrectedPairContext<'_, '_, R, V> {
    fn read_pair(&self) -> ReadPair<'_> {
        let pair = &self.corrected_pair;
        let fwd = crate::SequenceRead::from_views(
            pair.id(),
            pair.fwd_sequence(),
            pair.fwd_quality_scores(),
        );
        let rev = crate::SequenceRead::from_views(
            pair.id(),
            pair.rev_sequence(),
            pair.rev_quality_scores(),
        );
        ReadPair::from_views(fwd, rev)
    }
}

impl<R, V> HasPairOverlap for CorrectedPairContext<'_, '_, R, V> {
    fn pair_overlap(&self) -> Result<PairOverlap<'_>> {
        let overlap = match &self.overlap_outcome {
            OverlapOutcome::Found(overlap) => overlap,
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                return Err(OverlapError::NoOverlapFound.into());
            },
        };
        let bounds = overlap.bounds();
        let prepared = PreparedPair::from_read_pair(self.read_pair());

        PairOverlap::from_prepared(prepared, bounds)
    }
}

impl HasPairOverlap for ValidatedOverlap<'_> {
    fn pair_overlap(&self) -> Result<PairOverlap<'_>> {
        Ok(self.overlap().clone())
    }
}

impl HasConsensusRecord for MergedRead {
    fn consensus_id(&self) -> &str {
        self.id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.sequence()
    }

    fn consensus_quality_score_bytes(&self) -> &[u8] {
        self.quality_score_bytes()
    }
}

impl<V, C> HasConsensusRecord for MergeContext<'_, V, C> {
    fn consensus_id(&self) -> &str {
        self.merged_ref().id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.merged_ref().sequence()
    }

    fn consensus_quality_score_bytes(&self) -> &[u8] {
        self.merged_ref().quality_score_bytes()
    }
}

impl<V> HasConsensusRecord for CorrectedMergeContext<'_, V> {
    fn consensus_id(&self) -> &str {
        self.corrected_merged.id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.corrected_merged.sequence_bytes()
    }

    fn consensus_quality_score_bytes(&self) -> &[u8] {
        self.corrected_merged.quality_score_bytes()
    }
}

impl HasCorrectionWindow for MergedRead {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        Ok(CorrectionWindow::new(
            self.provenance().fwd_overlap_seq(),
            self.provenance().fwd_overlap_quality_score_bytes(),
            self.provenance().rev_overlap_seq(),
            self.provenance().rev_overlap_quality_score_bytes(),
        ))
    }
}

impl<V, C> HasCorrectionWindow for MergeContext<'_, V, C> {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        self.merged_ref().correction_window()
    }
}

impl<'asm, 'pair, R, V> HasCorrectionWindow
    for PairContext<
        'asm,
        'pair,
        R,
        super::typestate::HasOverlap,
        V,
        super::typestate::Unmerged,
        super::typestate::Uncorrected,
    >
where
    PairContext<
        'asm,
        'pair,
        R,
        super::typestate::HasOverlap,
        V,
        super::typestate::Unmerged,
        super::typestate::Uncorrected,
    >: HasReadPair + HasPairOverlap,
{
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        let overlap = match self.overlap_outcome() {
            OverlapOutcome::Found(overlap) => overlap,
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                return Err(OverlapError::NoOverlapFound.into());
            },
        };
        Ok(CorrectionWindow::from_overlap(overlap))
    }
}

impl HasConsensusRecord for CorrectedMergedRead {
    fn consensus_id(&self) -> &str {
        self.id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.sequence_bytes()
    }

    fn consensus_quality_score_bytes(&self) -> &[u8] {
        self.quality_score_bytes()
    }
}

impl<C> HasValidationMetrics for MergeContext<'_, super::typestate::Validated, C> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated merged contexts must retain validation metrics")
    }
}

impl HasValidationMetrics for CorrectedMergeContext<'_, super::typestate::Validated> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated corrected merged contexts must retain validation metrics")
    }
}

impl<R, M, C> HasValidationMetrics
    for PairContext<'_, '_, R, super::typestate::HasOverlap, super::typestate::Validated, M, C>
{
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated contexts must retain validation metrics")
    }
}

impl<R> HasValidationMetrics for CorrectedPairContext<'_, '_, R, super::typestate::Validated> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated corrected pair contexts must retain validation metrics")
    }
}

impl HasValidationMetrics for ValidatedOverlap<'_> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        ValidatedOverlap::validation_metrics(self)
    }
}
