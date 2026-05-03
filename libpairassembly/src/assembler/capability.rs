//! Internal capability traits for post-overlap operation contracts.

use crate::{
    OverlapParams, OverlapValidator, PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectedOrientedPair, CorrectionParams},
    errors::OverlapError,
    merge::{MergeParams, MergedRead},
    overlap::{HasOrientedPairSlices, OrientedPairSlices, OverlapBounds},
    validate::{ValidatedOverlap, ValidationMetrics},
};

use super::{
    Assembler,
    context::{CorrectedMergeContext, CorrectedPairContext, MergeContext, PairContext},
    typestate::{Corrected, Merged, OverlapFound, OverlapStateStorage, Unmerged},
};

pub(crate) mod private {
    pub(crate) trait Sealed {}
}

/// Internal marker trait for state/output carriers participating in assembler DAG operations.
pub(crate) trait PairState: private::Sealed {}

/// Capability for fluent assembler contexts that carry the full operation typestate tuple.
pub(crate) trait AssemblyContext: PairState {
    type OverlapState;
    type ValidationState;
    type MergeState;
    type CorrectionState;

    fn assembler(&self) -> &Assembler;

    #[inline]
    fn overlap_params(&self) -> &OverlapParams {
        self.assembler().overlap_params()
    }

    #[inline]
    fn validator(&self) -> &OverlapValidator {
        self.assembler().validator()
    }

    #[inline]
    fn correction_params(&self) -> CorrectionParams {
        self.assembler().correction_params()
    }

    #[inline]
    fn merge_params(&self) -> MergeParams {
        self.assembler().merge_params()
    }
}

/// Capability for exposing canonical oriented overlap slices for the current pair state.
pub(crate) trait HasPairOverlap: PairState {
    type Slices: HasOrientedPairSlices;

    fn pair_slices(&self) -> Result<&Self::Slices>;
    fn overlap_bounds(&self) -> Result<OverlapBounds>;

    fn validate_overlap_bounds(&self) -> Result<()> {
        self.pair_slices()?
            .validate_overlap_bounds(self.overlap_bounds()?)
    }

    fn overlap_windows(&self) -> Result<(&[u8], &[u8])> {
        let slices = self.pair_slices()?;
        let bounds = self.overlap_bounds()?;

        Ok((
            &slices.forward_sequence()[bounds.forward_range()],
            &slices.reverse_sequence_rc()[bounds.reverse_range()],
        ))
    }

    fn overlap_quality_windows(&self) -> Result<(&[u8], &[u8])> {
        let slices = self.pair_slices()?;
        let bounds = self.overlap_bounds()?;

        Ok((
            &slices.forward_quality_scores()[bounds.forward_range()],
            &slices.reverse_quality_scores_rc()[bounds.reverse_range()],
        ))
    }
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

impl<R, O, V, M, C> private::Sealed for PairContext<'_, '_, R, O, V, M, C> where
    O: for<'pair> OverlapStateStorage<'pair>
{
}
impl<R, O, V, M, C> PairState for PairContext<'_, '_, R, O, V, M, C> where
    O: for<'pair> OverlapStateStorage<'pair>
{
}

impl<R, O, V, M, C> AssemblyContext for PairContext<'_, '_, R, O, V, M, C>
where
    O: for<'pair> OverlapStateStorage<'pair>,
{
    type OverlapState = O;
    type ValidationState = V;
    type MergeState = M;
    type CorrectionState = C;

    fn assembler(&self) -> &Assembler {
        self.assembler
    }
}

impl private::Sealed for ValidatedOverlap<'_> {}
impl PairState for ValidatedOverlap<'_> {}

impl private::Sealed for PairOverlap<'_> {}
impl PairState for PairOverlap<'_> {}

impl private::Sealed for MergedRead {}
impl PairState for MergedRead {}

impl<V, C> private::Sealed for MergeContext<'_, '_, V, C> {}
impl<V, C> PairState for MergeContext<'_, '_, V, C> {}

impl<V, C> AssemblyContext for MergeContext<'_, '_, V, C> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Merged;
    type CorrectionState = C;

    fn assembler(&self) -> &Assembler {
        self.assembler
    }
}

impl<R, V> private::Sealed for CorrectedPairContext<'_, '_, R, V> {}
impl<R, V> PairState for CorrectedPairContext<'_, '_, R, V> {}

impl<R, V> AssemblyContext for CorrectedPairContext<'_, '_, R, V> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Unmerged;
    type CorrectionState = Corrected;

    fn assembler(&self) -> &Assembler {
        self.assembler
    }
}

impl<V> private::Sealed for CorrectedMergeContext<'_, V> {}
impl<V> PairState for CorrectedMergeContext<'_, V> {}

impl<V> AssemblyContext for CorrectedMergeContext<'_, V> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Merged;
    type CorrectionState = Corrected;

    fn assembler(&self) -> &Assembler {
        self.assembler
    }
}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl<'pair, R, V, M, C> HasPairOverlap for PairContext<'_, 'pair, R, OverlapFound, V, M, C> {
    type Slices = OrientedPairSlices<'pair>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.overlap().oriented_slices())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.overlap().bounds())
    }
}

impl<R, V> HasPairOverlap for CorrectedPairContext<'_, '_, R, V> {
    type Slices = CorrectedOrientedPair;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(&self.corrected_pair)
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.corrected_pair.overlap_bounds())
    }
}

impl<'a> HasPairOverlap for PairOverlap<'a> {
    type Slices = OrientedPairSlices<'a>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.oriented_slices())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.bounds())
    }
}

impl<'a> HasPairOverlap for ValidatedOverlap<'a> {
    type Slices = OrientedPairSlices<'a>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.overlap().oriented_slices())
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
    for PairContext<'_, '_, R, super::typestate::OverlapFound, super::typestate::Validated, M, C>
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
