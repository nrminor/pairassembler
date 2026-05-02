//! Transition legality markers and operation implementations.

use std::marker::PhantomData;

use crate::{
    OwnedSequenceRead, Result,
    correct::{CorrectedMergedRead, CorrectedPairEvidence},
    errors::OverlapError,
};

use super::{
    CorrectedContext, CorrectedMergeContext, CorrectedMergedContext, CorrectedPairContext,
    MergeContext, MergedContext, OverlapContext, PairContext, PairReady, SeqRecordView,
    ValidatedContext, ValidatedCorrectedContext, ValidatedCorrectedMergedContext,
    ValidatedMergedContext,
    capability::{HasMergeableOverlap, HasPairOverlap, HasValidationMetrics},
    context::OverlapOutcome,
    typestate::{Corrected, HasOverlap, NoOverlap, Uncorrected, Unmerged, Unvalidated, Validated},
};

#[derive(Debug, Clone)]
pub(crate) struct CanTuple<O, V, M, C>(PhantomData<(O, V, M, C)>);

pub(crate) trait CanOverlap {}
pub(crate) trait CanValidate {}
pub(crate) trait CanMerge {}
pub(crate) trait CanCorrect {}
impl CanOverlap for CanTuple<NoOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanValidate for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanValidate for CanTuple<HasOverlap, Unvalidated, Unmerged, Corrected> {}
impl CanMerge for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanMerge for CanTuple<HasOverlap, Validated, Unmerged, Uncorrected> {}
impl CanMerge for CanTuple<HasOverlap, Unvalidated, Unmerged, Corrected> {}
impl CanMerge for CanTuple<HasOverlap, Validated, Unmerged, Corrected> {}
impl CanCorrect for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanCorrect for CanTuple<HasOverlap, Validated, Unmerged, Uncorrected> {}

pub(crate) trait OverlapOp {
    type Out;
    fn overlap(self) -> Result<Self::Out>;
}

pub(crate) trait ValidateOp {
    type Out;
    fn validate(self) -> Result<Self::Out>;
}

pub(crate) trait ValidatePredicateOp {
    fn is_valid(&self) -> Result<bool>;
}

pub(crate) trait MergeOp {
    type Out;
    fn merge(self) -> Result<Self::Out>;
}

pub(crate) trait CorrectOp {
    type Out;
    fn correct(self) -> Result<Self::Out>;
}

impl<'asm, 'pair, R, O, V, M, C> OverlapOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanOverlap,
{
    type Out = OverlapContext<'asm, 'pair, R>;

    fn overlap(self) -> Result<Self::Out> {
        let PairContext {
            assembler,
            input,
            read_pair,
            ..
        } = self;
        let prepared = read_pair.prepare_for_overlap();
        let overlap_outcome =
            match prepared.scan_for_overlap_span_both(assembler.overlap_params())? {
                Some(overlap_span) => {
                    OverlapOutcome::Found(crate::PairOverlap::from_span(prepared, overlap_span)?)
                },
                None => OverlapOutcome::Missing,
            };

        Ok(PairContext {
            assembler,
            input,
            read_pair,
            overlap_outcome,
            validation_metrics: None,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R, O, V, M, C> ValidateOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanValidate,
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap,
{
    type Out = ValidatedContext<'asm, 'pair, R>;

    fn validate(self) -> Result<Self::Out> {
        let PairContext {
            assembler,
            input,
            read_pair,
            overlap_outcome,
            ..
        } = self;

        match overlap_outcome {
            OverlapOutcome::Found(overlap) => {
                let metrics = assembler.validator().assess(&overlap)?;
                Ok(PairContext {
                    assembler,
                    input,
                    read_pair,
                    overlap_outcome: OverlapOutcome::Found(overlap),
                    validation_metrics: Some(metrics),
                    _marker: PhantomData,
                })
            },
            OverlapOutcome::Missing | OverlapOutcome::Unknown => Ok(PairContext {
                assembler,
                input,
                read_pair,
                overlap_outcome: OverlapOutcome::Missing,
                validation_metrics: None,
                _marker: PhantomData,
            }),
        }
    }
}

impl<'asm, 'pair, R> ValidateOp for CorrectedPairContext<'asm, 'pair, R, Unvalidated>
where
    R: SeqRecordView,
{
    type Out = ValidatedCorrectedContext<'asm, 'pair, R>;

    fn validate(self) -> Result<Self::Out> {
        let metrics = {
            let validator = self.validator();

            validator.assess(&self)?
        };

        let CorrectedPairContext {
            assembler,
            input,
            corrected_pair,
            overlap_bounds,
            validation_metrics: _,
            _marker: _,
        } = self;

        Ok(CorrectedPairContext {
            assembler,
            input,
            corrected_pair,
            overlap_bounds,
            validation_metrics: Some(metrics),
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R, O, V, M, C> ValidatePredicateOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap,
{
    fn is_valid(&self) -> Result<bool> {
        if self.validation_metrics_ref().is_some() {
            return Ok(true);
        }

        match self.overlap_outcome() {
            OverlapOutcome::Found(_) => Ok(self.validator().assess(self).is_ok()),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => Ok(false),
        }
    }
}

impl<'asm, 'pair, R, V> MergeOp
    for PairContext<'asm, 'pair, R, HasOverlap, V, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, V, Unmerged, Uncorrected>: CanMerge,
    Self: HasMergeableOverlap,
{
    type Out = MergeContext<'asm, 'pair, V, Uncorrected>;

    fn merge(self) -> Result<Self::Out> {
        let consensus = self.merge_consensus()?;

        let PairContext {
            assembler,
            overlap_outcome,
            validation_metrics,
            ..
        } = self;

        match overlap_outcome {
            OverlapOutcome::Found(overlap) => Ok(MergeContext {
                assembler,
                consensus,
                overlap,
                validation_metrics,
                _marker: PhantomData,
            }),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }
}

impl<'asm, 'pair, R, V> MergeOp for CorrectedPairContext<'asm, 'pair, R, V>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, V, Unmerged, Corrected>: CanMerge,
{
    type Out = CorrectedMergeContext<'asm, V>;

    fn merge(self) -> Result<Self::Out> {
        let CorrectedPairContext {
            assembler,
            corrected_pair,
            overlap_bounds,
            validation_metrics,
            ..
        } = self;
        let corrected_merged = {
            let consensus = corrected_pair.into_merged_consensus(overlap_bounds)?;
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
    for PairContext<'asm, 'pair, R, HasOverlap, V, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, V, Unmerged, Uncorrected>: CanCorrect,
    PairContext<'asm, 'pair, R, HasOverlap, V, Unmerged, Uncorrected>: HasPairOverlap,
{
    type Out = CorrectedPairContext<'asm, 'pair, R, Unvalidated>;

    fn correct(self) -> Result<Self::Out> {
        let correction = self.correction_params();
        let PairContext {
            assembler,
            input,
            read_pair,
            overlap_outcome,
            ..
        } = self;

        match overlap_outcome {
            OverlapOutcome::Found(overlap) => {
                let corrected_pair = CorrectedPairEvidence::correct_from_overlap_with(
                    &read_pair, &overlap, correction,
                );
                let overlap_bounds = overlap.bounds();
                Ok(super::CorrectedPairContext {
                    assembler,
                    input,
                    corrected_pair,
                    overlap_bounds,
                    validation_metrics: None,
                    _marker: PhantomData,
                })
            },
            OverlapOutcome::Missing | OverlapOutcome::Unknown => {
                Err(OverlapError::NoOverlapFound.into())
            },
        }
    }
}

impl<'asm, 'pair, V> CorrectOp for MergeContext<'asm, 'pair, V, Uncorrected> {
    type Out = CorrectedMergeContext<'asm, V>;

    fn correct(self) -> Result<Self::Out> {
        let correction = self.correction_params();
        let MergeContext {
            assembler,
            consensus,
            overlap,
            validation_metrics,
            ..
        } = self;
        let corrected_merged =
            CorrectedMergedRead::correct_consensus_with(consensus, &overlap, correction)?;
        Ok(super::CorrectedMergeContext {
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
    pub fn overlap(self) -> Result<OverlapContext<'asm, 'pair, R>> {
        OverlapOp::overlap(self)
    }

    /// Process this pair to completion using the parent assembler.
    ///
    /// # Errors
    ///
    /// Returns any pipeline error encountered while processing this pair.
    pub fn process(self) -> Result<OwnedSequenceRead> {
        let corrected =
            MergeOp::merge(ValidateOp::validate(OverlapOp::overlap(self)?)?)?.correct()?;
        corrected.into_owned_read()
    }
}

impl<'asm, 'pair, V> MergeContext<'asm, 'pair, V, Uncorrected> {
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

    /// Merge this overlap using the current unvalidated pair evidence.
    ///
    /// # Errors
    ///
    /// Returns an error if merge fails.
    pub fn merge(self) -> Result<MergedContext<'asm, 'pair>> {
        MergeOp::merge(self)
    }

    /// Correct both mates directly from overlap evidence.
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
    /// Validate corrected pair evidence with configured validator and enter validated corrected context.
    ///
    /// # Errors
    ///
    /// Returns an error if corrected overlap validation fails.
    pub fn validate(self) -> Result<ValidatedCorrectedContext<'asm, 'pair, R>> {
        ValidateOp::validate(self)
    }

    /// Merge corrected pair evidence in its current unvalidated state.
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
    /// Merge corrected pair evidence after corrected-evidence validation.
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

    /// Correct both mates directly from overlap evidence.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap-derived correction fails.
    pub fn correct(self) -> Result<CorrectedContext<'asm, 'pair, R>> {
        CorrectOp::correct(self)
    }
}
