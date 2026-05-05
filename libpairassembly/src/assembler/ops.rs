//! Transition legality markers and operation implementations.

use std::marker::PhantomData;

use crate::{
    PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectionParams, OverlapCorrector},
    merge::OverlapMerger,
    overlap::{OrientedPairSlices, OverlapFinder},
    read::OwnedSequenceRead,
    validate::ValidationMetrics,
};

use super::{
    CorrectedContext, CorrectedMergedContext, MergedContext, NoOverlapContext, OverlapContext,
    OverlapOutcome, OverlapSearch, PairReady, SeqRecordView, ValidatedContext,
    ValidatedCorrectedContext, ValidatedCorrectedMergedContext, ValidatedMergedContext,
    capability::{AssemblyContext, HasPairOverlap, HasValidationMetrics},
    context::{CorrectedMergeContext, CorrectedPairContext, MergeContext, PairContext},
    typestate::{
        Corrected, Merged, NoOverlapFound, OverlapFound, Uncorrected, Unmerged, Unvalidated,
        Validated,
    },
};

pub(crate) trait OverlapOp<'pair, 'scratch>: AssemblyContext + Sized {
    type Found: AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Uncorrected,
        > + HasPairOverlap;
    type NoOverlap: AssemblyContext<
            OverlapState = NoOverlapFound,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Uncorrected,
        >;

    fn oriented_slices(&self) -> OrientedPairSlices<'pair, 'scratch>;
    fn found_overlap(self, overlap: PairOverlap<'pair, 'scratch>) -> Self::Found;
    fn no_overlap_found(self) -> Self::NoOverlap;

    fn overlap(self) -> Result<OverlapOutcome<Self::Found, Self::NoOverlap>> {
        match OverlapFinder::new(self.overlap_params()).find_in_slices(self.oriented_slices())? {
            Some(overlap) => Ok(OverlapOutcome::Found(self.found_overlap(overlap))),
            None => Ok(OverlapOutcome::NoOverlap(self.no_overlap_found())),
        }
    }
}

pub(crate) trait ValidateOp:
    AssemblyContext<OverlapState = OverlapFound> + HasPairOverlap + Sized
{
    type Out: AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = Validated,
            MergeState = <Self as AssemblyContext>::MergeState,
            CorrectionState = <Self as AssemblyContext>::CorrectionState,
        > + HasPairOverlap
        + HasValidationMetrics;

    fn into_validated(self, metrics: ValidationMetrics) -> Self::Out;

    fn validate(self) -> Result<Self::Out> {
        let metrics = self.validator().assess(&self)?;
        Ok(self.into_validated(metrics))
    }
}

pub(crate) trait ValidatePredicateOp:
    AssemblyContext<OverlapState = OverlapFound> + HasPairOverlap
{
    fn is_valid(&self) -> Result<bool>;
}

pub(crate) trait MergeOp:
    AssemblyContext<OverlapState = OverlapFound> + HasPairOverlap + Sized
{
    type Out: AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = <Self as AssemblyContext>::ValidationState,
            MergeState = Merged,
            CorrectionState = <Self as AssemblyContext>::CorrectionState,
        > + HasPairOverlap;

    fn merge(self) -> Result<Self::Out>;
}

pub(crate) trait CorrectOp: AssemblyContext<OverlapState = OverlapFound> + Sized {
    type Out: AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = Unvalidated,
            MergeState = <Self as AssemblyContext>::MergeState,
            CorrectionState = Corrected,
        >;

    fn correct_with_params(self, params: CorrectionParams) -> Result<Self::Out>;

    fn correct(self) -> Result<Self::Out> {
        let params = self.correction_params();
        self.correct_with_params(params)
    }
}

impl<'asm, 'pair, R> OverlapOp<'pair, 'asm> for PairReady<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    type Found = OverlapContext<'asm, 'pair, R>;
    type NoOverlap = NoOverlapContext<'asm, 'pair, R>;

    fn oriented_slices(&self) -> OrientedPairSlices<'pair, 'asm> {
        self.overlap
    }

    fn found_overlap(self, overlap: PairOverlap<'pair, 'asm>) -> Self::Found {
        let PairContext {
            config,
            input,
            read_pair,
            ..
        } = self;

        PairContext {
            config,
            input,
            read_pair,
            overlap,
            validation_metrics: None,
            _marker: PhantomData,
        }
    }

    fn no_overlap_found(self) -> Self::NoOverlap {
        let PairContext {
            config,
            input,
            read_pair,
            ..
        } = self;

        PairContext {
            config,
            input,
            read_pair,
            overlap: (),
            validation_metrics: None,
            _marker: PhantomData,
        }
    }
}

impl<'asm, 'pair, R, V> ValidateOp
    for PairContext<'asm, 'pair, R, OverlapFound, V, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    PairContext<'asm, 'pair, R, OverlapFound, V, Unmerged, Uncorrected>: HasPairOverlap,
{
    type Out = ValidatedContext<'asm, 'pair, R>;

    fn into_validated(self, metrics: ValidationMetrics) -> Self::Out {
        let PairContext {
            config,
            input,
            read_pair,
            overlap,
            ..
        } = self;

        PairContext {
            config,
            input,
            read_pair,
            overlap,
            validation_metrics: Some(metrics),
            _marker: PhantomData,
        }
    }
}

impl<'asm, 'pair, R> ValidateOp for CorrectedPairContext<'asm, 'pair, R, Unvalidated>
where
    R: SeqRecordView,
{
    type Out = ValidatedCorrectedContext<'asm, 'pair, R>;

    fn into_validated(self, metrics: ValidationMetrics) -> Self::Out {
        let CorrectedPairContext {
            config,
            input,
            corrected_pair,
            validation_metrics: _,
            _marker: _,
        } = self;

        CorrectedPairContext {
            config,
            input,
            corrected_pair,
            validation_metrics: Some(metrics),
            _marker: PhantomData,
        }
    }
}

impl<'asm, 'pair, C> ValidateOp for MergeContext<'asm, 'pair, Unvalidated, C> {
    type Out = MergeContext<'asm, 'pair, Validated, C>;

    fn into_validated(self, metrics: ValidationMetrics) -> Self::Out {
        let MergeContext {
            config,
            consensus,
            overlap,
            ..
        } = self;

        MergeContext {
            config,
            consensus,
            overlap,
            validation_metrics: Some(metrics),
            _marker: PhantomData,
        }
    }
}

impl<'asm, 'pair> ValidateOp for CorrectedMergeContext<'asm, 'pair, Unvalidated> {
    type Out = CorrectedMergeContext<'asm, 'pair, Validated>;

    fn into_validated(self, metrics: ValidationMetrics) -> Self::Out {
        let CorrectedMergeContext {
            config,
            corrected_merged,
            corrected_pair,
            ..
        } = self;

        CorrectedMergeContext {
            config,
            corrected_merged,
            corrected_pair,
            validation_metrics: Some(metrics),
            _marker: PhantomData,
        }
    }
}

impl<'asm, 'pair, R, V, M, C> ValidatePredicateOp
    for PairContext<'asm, 'pair, R, OverlapFound, V, M, C>
where
    PairContext<'asm, 'pair, R, OverlapFound, V, M, C>: HasPairOverlap,
{
    fn is_valid(&self) -> Result<bool> {
        if self.validation_metrics_ref().is_some() {
            return Ok(true);
        }

        Ok(self.validator().assess(self).is_ok())
    }
}

impl<'asm, 'pair, R, V> MergeOp
    for PairContext<'asm, 'pair, R, OverlapFound, V, Unmerged, Uncorrected>
{
    type Out = MergeContext<'asm, 'pair, V, Uncorrected>;

    fn merge(self) -> Result<Self::Out> {
        let consensus = OverlapMerger::new(self.merge_params()).merge_consensus(&self)?;

        let PairContext {
            config,
            overlap,
            validation_metrics,
            ..
        } = self;

        Ok(MergeContext {
            config,
            consensus,
            overlap,
            validation_metrics,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R, V> MergeOp for CorrectedPairContext<'asm, 'pair, R, V> {
    type Out = CorrectedMergeContext<'asm, 'pair, V>;

    fn merge(self) -> Result<Self::Out> {
        let CorrectedPairContext {
            config,
            corrected_pair,
            validation_metrics,
            ..
        } = self;
        let corrected_merged =
            CorrectedMergedRead::try_from(corrected_pair.to_merged_consensus()?)?;

        Ok(CorrectedMergeContext {
            config,
            corrected_merged,
            corrected_pair,
            validation_metrics,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R, V> CorrectOp
    for PairContext<'asm, 'pair, R, OverlapFound, V, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    PairContext<'asm, 'pair, R, OverlapFound, V, Unmerged, Uncorrected>: HasPairOverlap,
{
    type Out = CorrectedPairContext<'asm, 'pair, R, Unvalidated>;

    fn correct_with_params(self, correction: CorrectionParams) -> Result<Self::Out> {
        let corrected_pair = OverlapCorrector::new(correction).correct_pair_overlap(&self)?;

        let PairContext { config, input, .. } = self;

        Ok(CorrectedPairContext {
            config,
            input,
            corrected_pair,
            validation_metrics: None,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, V> CorrectOp for MergeContext<'asm, 'pair, V, Uncorrected> {
    type Out = CorrectedMergeContext<'asm, 'pair, Unvalidated>;

    fn correct_with_params(self, correction: CorrectionParams) -> Result<Self::Out> {
        let corrector = OverlapCorrector::new(correction);
        let corrected_pair = corrector.correct_pair_overlap(&self)?;

        let MergeContext {
            config,
            consensus,
            overlap,
            ..
        } = self;
        let corrected_merged = corrector.correct_merged_consensus(consensus, &overlap)?;
        Ok(CorrectedMergeContext {
            config,
            corrected_merged,
            corrected_pair,
            validation_metrics: None,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R> PairReady<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    /// Detect overlap and enter overlap context.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap discovery fails.
    pub fn find_overlap(self) -> Result<OverlapSearch<'asm, 'pair, R>> {
        OverlapOp::overlap(self)
    }

    /// Process this pair to completion using the parent assembler.
    ///
    /// # Errors
    ///
    /// Returns any pipeline error encountered while processing this pair.
    pub fn process(self) -> Result<Option<OwnedSequenceRead>> {
        OverlapOp::overlap(self)?.and_then_found(|overlap| {
            let corrected = MergeOp::merge(ValidateOp::validate(overlap)?)?.correct()?;
            corrected.into_owned_read()
        })
    }
}

impl<'asm, 'pair, V> MergeContext<'asm, 'pair, V, Uncorrected> {
    /// Correct this merged artifact using the configured correction policy.
    ///
    /// # Errors
    ///
    /// Returns an error if correction fails for the merged artifact.
    pub fn correct(self) -> Result<CorrectedMergeContext<'asm, 'pair, Unvalidated>> {
        CorrectOp::correct(self)
    }
}

impl<'asm, 'pair, C> MergeContext<'asm, 'pair, Unvalidated, C> {
    /// Validate the overlap evidence retained by this merged context.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation fails.
    pub fn validate(self) -> Result<MergeContext<'asm, 'pair, Validated, C>> {
        ValidateOp::validate(self)
    }
}

impl<'asm, 'pair> CorrectedMergeContext<'asm, 'pair, Unvalidated> {
    /// Validate the corrected overlap evidence retained by this corrected merged context.
    ///
    /// # Errors
    ///
    /// Returns an error if corrected-overlap validation fails.
    pub fn validate(self) -> Result<CorrectedMergeContext<'asm, 'pair, Validated>> {
        ValidateOp::validate(self)
    }
}

impl<'asm, 'pair, R> OverlapContext<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    /// Check overlap validity without transitioning typestate.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation cannot be evaluated.
    pub fn is_valid(&self) -> Result<bool> {
        ValidatePredicateOp::is_valid(self)
    }

    /// Validate overlap with configured validator and enter validated context.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation fails.
    pub fn validate(self) -> Result<ValidatedContext<'asm, 'pair, R>> {
        ValidateOp::validate(self)
    }

    /// Merge this overlap using the current unvalidated pair slices.
    ///
    /// # Errors
    ///
    /// Returns an error if merge fails.
    pub fn merge(self) -> Result<MergedContext<'asm, 'pair>> {
        MergeOp::merge(self)
    }

    /// Correct both mates directly from overlap slices.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap-derived correction fails.
    pub fn correct(self) -> Result<CorrectedContext<'asm, 'pair, R>> {
        CorrectOp::correct(self)
    }
}

impl<'asm, 'pair, R> CorrectedContext<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    /// Validate corrected pair slices with configured validator and enter validated corrected context.
    ///
    /// # Errors
    ///
    /// Returns an error if corrected overlap validation fails.
    pub fn validate(self) -> Result<ValidatedCorrectedContext<'asm, 'pair, R>> {
        ValidateOp::validate(self)
    }

    /// Merge corrected pair slices in its current unvalidated state.
    ///
    /// # Errors
    ///
    /// Returns an error if merge projection or consensus construction fails.
    pub fn merge(self) -> Result<CorrectedMergedContext<'asm, 'pair>> {
        MergeOp::merge(self)
    }
}

impl<'asm, 'pair, R> ValidatedCorrectedContext<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    /// Merge corrected pair slices after corrected-slice validation.
    ///
    /// # Errors
    ///
    /// Returns an error if merge projection or consensus construction fails.
    pub fn merge(self) -> Result<ValidatedCorrectedMergedContext<'asm, 'pair>> {
        MergeOp::merge(self)
    }
}

impl<'asm, 'pair, R> ValidatedContext<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    /// Check overlap validity without transitioning typestate.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation cannot be evaluated.
    pub fn is_valid(&self) -> Result<bool> {
        ValidatePredicateOp::is_valid(self)
    }

    /// Merge this pair using the checked path.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation or merge fails.
    pub fn merge(self) -> Result<ValidatedMergedContext<'asm, 'pair>> {
        MergeOp::merge(self)
    }

    /// Correct both mates directly from overlap slices.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap-derived correction fails.
    pub fn correct(self) -> Result<CorrectedContext<'asm, 'pair, R>> {
        CorrectOp::correct(self)
    }
}
