//! Internal capability traits for post-overlap operation contracts.

use crate::{
    PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectedPairEvidence, CorrectionWindow},
    errors::OverlapError,
    merge::{MergeView, MergedConsensus, MergedRead},
    overlap::{HasOrientedPairEvidence, OverlapBounds, PreparedPair},
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
    type Evidence: HasOrientedPairEvidence;

    fn pair_evidence(&self) -> Result<&Self::Evidence>;
    fn overlap_bounds(&self) -> Result<OverlapBounds>;
}

/// Capability for exposing normalized merge-ready overlap views.
pub(crate) trait HasMergeableOverlap: HasPairOverlap {
    fn merge_view(&self) -> Result<MergeView<'_>> {
        MergeView::from_oriented_evidence(self.pair_evidence()?, self.overlap_bounds()?)
    }

    fn merge_consensus(&self) -> Result<MergedConsensus> {
        MergedConsensus::try_from_merge_view(self.merge_view()?)
    }
}

impl<T> HasMergeableOverlap for T where T: HasPairOverlap {}

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

impl private::Sealed for PairOverlap<'_> {}
impl PairState for PairOverlap<'_> {}

impl private::Sealed for MergedRead {}
impl PairState for MergedRead {}

impl<V, C> private::Sealed for MergeContext<'_, '_, V, C> {}
impl<V, C> PairState for MergeContext<'_, '_, V, C> {}

impl<R, V> private::Sealed for CorrectedPairContext<'_, '_, R, V> {}
impl<R, V> PairState for CorrectedPairContext<'_, '_, R, V> {}

impl<V> private::Sealed for CorrectedMergeContext<'_, V> {}
impl<V> PairState for CorrectedMergeContext<'_, V> {}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl<'pair, R, O, V, M, C> HasPairOverlap for PairContext<'_, 'pair, R, O, V, M, C> {
    type Evidence = PreparedPair<'pair>;

    fn pair_evidence(&self) -> Result<&Self::Evidence> {
        match self.overlap_outcome() {
            OverlapOutcome::Found(overlap) => Ok(overlap.prepared_evidence()),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        match self.overlap_outcome() {
            OverlapOutcome::Found(overlap) => Ok(overlap.bounds()),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }
}

impl<R, V> HasPairOverlap for CorrectedPairContext<'_, '_, R, V> {
    type Evidence = CorrectedPairEvidence;

    fn pair_evidence(&self) -> Result<&Self::Evidence> {
        Ok(&self.corrected_pair)
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.overlap_bounds)
    }
}

impl<'a> HasPairOverlap for PairOverlap<'a> {
    type Evidence = PreparedPair<'a>;

    fn pair_evidence(&self) -> Result<&Self::Evidence> {
        Ok(self.prepared_evidence())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.bounds())
    }
}

impl<'a> HasPairOverlap for ValidatedOverlap<'a> {
    type Evidence = PreparedPair<'a>;

    fn pair_evidence(&self) -> Result<&Self::Evidence> {
        Ok(self.overlap().prepared_evidence())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.overlap().bounds())
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

impl<V, C> HasConsensusRecord for MergeContext<'_, '_, V, C> {
    fn consensus_id(&self) -> &str {
        self.consensus_ref().id()
    }

    fn consensus_seq(&self) -> &[u8] {
        self.consensus_ref().sequence()
    }

    fn consensus_quality_score_bytes(&self) -> &[u8] {
        self.consensus_ref().quality_score_bytes()
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

impl<V, C> HasCorrectionWindow for MergeContext<'_, '_, V, C> {
    fn correction_window(&self) -> Result<CorrectionWindow<'_>> {
        Ok(CorrectionWindow::from_overlap(self.overlap_ref()))
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
    >: HasPairOverlap,
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

impl<C> HasValidationMetrics for MergeContext<'_, '_, super::typestate::Validated, C> {
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
