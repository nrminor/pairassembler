//! Transition legality markers and operation implementations.

use std::marker::PhantomData;

use crate::{
    PairOverlap, Result,
    correct::{CorrectedMergedRead, CorrectionParams, OverlapCorrector},
    merge::OverlapMerger,
    overlap::OverlapFinder,
    read::{OwnedSequenceRead, ReadPair},
    validate::ValidationMetrics,
};

use super::{
    CorrectedContext, CorrectedMergedContext, MergedContext, NoOverlapContext, OverlapContext,
    OverlapOutcome, OverlapSearch, PairReady, SeqRecordView, ValidatedContext,
    ValidatedCorrectedContext, ValidatedCorrectedMergedContext, ValidatedMergedContext,
    capability::{AssemblyContext, HasConsensusRecord, HasPairOverlap, HasValidationMetrics},
    context::{CorrectedMergeContext, CorrectedPairContext, MergeContext, PairContext},
    typestate::{
        Corrected, Merged, NoOverlapFound, OverlapFound, Uncorrected, Unmerged, Unvalidated,
        Validated,
    },
};

pub(crate) trait OverlapOp<'pair>: AssemblyContext + Sized {
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

    fn read_pair(&self) -> ReadPair<'pair>;
    fn found_overlap(self, overlap: PairOverlap<'pair>) -> Self::Found;
    fn no_overlap_found(self) -> Self::NoOverlap;

    fn overlap(self) -> Result<OverlapOutcome<Self::Found, Self::NoOverlap>> {
        match OverlapFinder::new(self.overlap_params()).find(self.read_pair())? {
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
        > + HasConsensusRecord;

    fn merge(self) -> Result<Self::Out>;
}

pub(crate) trait CorrectOp: AssemblyContext<OverlapState = OverlapFound> + Sized {
    type Out: AssemblyContext<
            OverlapState = OverlapFound,
            MergeState = <Self as AssemblyContext>::MergeState,
            CorrectionState = Corrected,
        >;

    fn correct_with_params(self, params: CorrectionParams) -> Result<Self::Out>;

    fn correct(self) -> Result<Self::Out> {
        let params = self.correction_params();
        self.correct_with_params(params)
    }
}

impl<'asm, 'pair, R> OverlapOp<'pair> for PairReady<'asm, 'pair, R>
where
    R: SeqRecordView,
{
    type Found = OverlapContext<'asm, 'pair, R>;
    type NoOverlap = NoOverlapContext<'asm, 'pair, R>;

    fn read_pair(&self) -> ReadPair<'pair> {
        self.read_pair
    }

    fn found_overlap(self, overlap: PairOverlap<'pair>) -> Self::Found {
        let PairContext {
            assembler,
            input,
            read_pair,
            ..
        } = self;

        PairContext {
            assembler,
            input,
            read_pair,
            overlap,
            validation_metrics: None,
            _marker: PhantomData,
        }
    }

    fn no_overlap_found(self) -> Self::NoOverlap {
        let PairContext {
            assembler,
            input,
            read_pair,
            ..
        } = self;

        PairContext {
            assembler,
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
            assembler,
            input,
            read_pair,
            overlap,
            ..
        } = self;

        PairContext {
            assembler,
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
            assembler,
            input,
            corrected_pair,
            validation_metrics: _,
            _marker: _,
        } = self;

        CorrectedPairContext {
            assembler,
            input,
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
            assembler,
            overlap,
            validation_metrics,
            ..
        } = self;

        Ok(MergeContext {
            assembler,
            consensus,
            overlap,
            validation_metrics,
            _marker: PhantomData,
        })
    }
}

impl<'asm, R, V> MergeOp for CorrectedPairContext<'asm, '_, R, V> {
    type Out = CorrectedMergeContext<'asm, V>;

    fn merge(self) -> Result<Self::Out> {
        let CorrectedPairContext {
            assembler,
            corrected_pair,
            validation_metrics,
            ..
        } = self;
        let corrected_merged = {
            let consensus = corrected_pair.into_merged_consensus()?;
            CorrectedMergedRead::try_from(consensus)?
        };

        Ok(CorrectedMergeContext {
            assembler,
            corrected_merged,
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

        let PairContext {
            assembler, input, ..
        } = self;

        Ok(CorrectedPairContext {
            assembler,
            input,
            corrected_pair,
            validation_metrics: None,
            _marker: PhantomData,
        })
    }
}

impl<'asm, V> CorrectOp for MergeContext<'asm, '_, V, Uncorrected> {
    type Out = CorrectedMergeContext<'asm, V>;

    fn correct_with_params(self, correction: CorrectionParams) -> Result<Self::Out> {
        let MergeContext {
            assembler,
            consensus,
            overlap,
            validation_metrics,
            ..
        } = self;
        let corrected_merged =
            OverlapCorrector::new(correction).correct_merged_consensus(consensus, &overlap)?;
        Ok(CorrectedMergeContext {
            assembler,
            corrected_merged,
            validation_metrics,
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

impl<'asm, V> MergeContext<'asm, '_, V, Uncorrected> {
    /// Correct this merged artifact using the configured correction policy.
    ///
    /// # Errors
    ///
    /// Returns an error if correction fails for the merged artifact.
    pub fn correct(self) -> Result<CorrectedMergeContext<'asm, V>> {
        CorrectOp::correct(self)
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
    pub fn merge(self) -> Result<CorrectedMergedContext<'asm>> {
        MergeOp::merge(self)
    }
}

impl<'asm, R> ValidatedCorrectedContext<'asm, '_, R>
where
    R: SeqRecordView,
{
    /// Merge corrected pair slices after corrected-slice validation.
    ///
    /// # Errors
    ///
    /// Returns an error if merge projection or consensus construction fails.
    pub fn merge(self) -> Result<ValidatedCorrectedMergedContext<'asm>> {
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
