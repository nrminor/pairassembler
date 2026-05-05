//! Internal capability traits for post-overlap operation contracts.

use crate::{
    OverlapParams, OverlapValidator, PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectedOrientedPair, CorrectionParams},
    merge::MergeParams,
    overlap::{HasOrientedPairSlices, OrientedPairSlices, OverlapBounds},
    validate::{ValidatedOverlap, ValidationMetrics},
};

use super::{
    AssemblerConfig,
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

    fn config(&self) -> &AssemblerConfig;

    #[inline]
    fn overlap_params(&self) -> &OverlapParams {
        self.config().overlap_params()
    }

    #[inline]
    fn validator(&self) -> &OverlapValidator {
        self.config().validator()
    }

    #[inline]
    fn correction_params(&self) -> CorrectionParams {
        self.config().correction_params()
    }

    #[inline]
    fn merge_params(&self) -> MergeParams {
        self.config().merge_params()
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
        self.validate_overlap_bounds()?;
        let slices = self.pair_slices()?;
        let bounds = self.overlap_bounds()?;

        Ok((
            &slices.forward_sequence()[bounds.forward_range()],
            &slices.reverse_sequence_rc()[bounds.reverse_range()],
        ))
    }

    fn overlap_quality_windows(&self) -> Result<(&[u8], &[u8])> {
        self.validate_overlap_bounds()?;
        let slices = self.pair_slices()?;
        let bounds = self.overlap_bounds()?;

        Ok((
            &slices.forward_quality_score_bytes()[bounds.forward_range()],
            &slices.reverse_quality_score_bytes_rc()[bounds.reverse_range()],
        ))
    }
}

/// Capability for exposing retained validation-stage metrics.
pub trait HasValidationMetrics {
    fn validation_metrics(&self) -> &ValidationMetrics;
}

impl<R, O, V, M, C> private::Sealed for PairContext<'_, '_, R, O, V, M, C> where
    O: for<'pair, 'scratch> OverlapStateStorage<'pair, 'scratch>
{
}
impl<R, O, V, M, C> PairState for PairContext<'_, '_, R, O, V, M, C> where
    O: for<'pair, 'scratch> OverlapStateStorage<'pair, 'scratch>
{
}

impl<R, O, V, M, C> AssemblyContext for PairContext<'_, '_, R, O, V, M, C>
where
    O: for<'pair, 'scratch> OverlapStateStorage<'pair, 'scratch>,
{
    type OverlapState = O;
    type ValidationState = V;
    type MergeState = M;
    type CorrectionState = C;

    fn config(&self) -> &AssemblerConfig {
        self.config
    }
}

impl private::Sealed for ValidatedOverlap<'_, '_> {}
impl PairState for ValidatedOverlap<'_, '_> {}

impl private::Sealed for PairOverlap<'_, '_> {}
impl PairState for PairOverlap<'_, '_> {}

impl<V, C> private::Sealed for MergeContext<'_, '_, V, C> {}
impl<V, C> PairState for MergeContext<'_, '_, V, C> {}

impl<V, C> AssemblyContext for MergeContext<'_, '_, V, C> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Merged;
    type CorrectionState = C;

    fn config(&self) -> &AssemblerConfig {
        self.config
    }
}

impl<R, V> private::Sealed for CorrectedPairContext<'_, '_, R, V> {}
impl<R, V> PairState for CorrectedPairContext<'_, '_, R, V> {}

impl<R, V> AssemblyContext for CorrectedPairContext<'_, '_, R, V> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Unmerged;
    type CorrectionState = Corrected;

    fn config(&self) -> &AssemblerConfig {
        self.config
    }
}

impl<V> private::Sealed for CorrectedMergeContext<'_, '_, V> {}
impl<V> PairState for CorrectedMergeContext<'_, '_, V> {}

impl<V> AssemblyContext for CorrectedMergeContext<'_, '_, V> {
    type OverlapState = OverlapFound;
    type ValidationState = V;
    type MergeState = Merged;
    type CorrectionState = Corrected;

    fn config(&self) -> &AssemblerConfig {
        self.config
    }
}

impl private::Sealed for CorrectedMergedRead {}
impl PairState for CorrectedMergedRead {}

impl<'scratch, 'pair, R, V, M, C> HasPairOverlap
    for PairContext<'scratch, 'pair, R, OverlapFound, V, M, C>
{
    type Slices = OrientedPairSlices<'pair, 'scratch>;

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

impl<'scratch, 'pair, V, C> HasPairOverlap for MergeContext<'scratch, 'pair, V, C> {
    type Slices = OrientedPairSlices<'pair, 'scratch>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.overlap.oriented_slices())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.overlap.bounds())
    }
}

impl<V> HasPairOverlap for CorrectedMergeContext<'_, '_, V> {
    type Slices = CorrectedOrientedPair;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(&self.corrected_pair)
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.corrected_pair.overlap_bounds())
    }
}

impl<'pair, 'scratch> HasPairOverlap for PairOverlap<'pair, 'scratch> {
    type Slices = OrientedPairSlices<'pair, 'scratch>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.oriented_slices())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.bounds())
    }
}

impl<'pair, 'scratch> HasPairOverlap for ValidatedOverlap<'pair, 'scratch> {
    type Slices = OrientedPairSlices<'pair, 'scratch>;

    fn pair_slices(&self) -> Result<&Self::Slices> {
        Ok(self.overlap().oriented_slices())
    }

    fn overlap_bounds(&self) -> Result<OverlapBounds> {
        Ok(self.overlap().bounds())
    }
}

impl<C> HasValidationMetrics for MergeContext<'_, '_, super::typestate::Validated, C> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        self.validation_metrics_ref()
            .expect("validated merged contexts must retain validation metrics")
    }
}

impl HasValidationMetrics for CorrectedMergeContext<'_, '_, super::typestate::Validated> {
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

impl HasValidationMetrics for ValidatedOverlap<'_, '_> {
    fn validation_metrics(&self) -> &ValidationMetrics {
        ValidatedOverlap::validation_metrics(self)
    }
}
