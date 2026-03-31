//! Transition legality markers and operation implementations.

use std::marker::PhantomData;

use crate::{
    Result,
    correct::{CorrectedMergedRead, CorrectedReadPair},
    errors::OverlapError,
    merge::{MergedRead, merge_from},
};

use super::{
    OverlapContext, PairContext, PairReady, SeqRecordView, ValidatedContext,
    capability::{HasPairOverlap, HasReadPair},
    context::{OverlapOutcome, OverlapSnapshot},
    typestate::{HasOverlap, NoOverlap, Uncorrected, Unmerged, Unvalidated, Validated},
};

#[derive(Debug, Clone)]
pub(crate) struct CanTuple<O, V, M, C>(PhantomData<(O, V, M, C)>);

pub(crate) trait CanOverlap {}
pub(crate) trait CanValidate {}
pub(crate) trait CanMerge {}
pub(crate) trait CanCorrectPair {}
pub(crate) trait CanCorrectPairUnchecked {}
pub(crate) trait CanCorrectMerged {
    fn into_corrected_merged(self) -> Result<CorrectedMergedRead>;
}

impl CanOverlap for CanTuple<NoOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanValidate for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanMerge for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanMerge for CanTuple<HasOverlap, Validated, Unmerged, Uncorrected> {}
impl CanCorrectPair for CanTuple<HasOverlap, Validated, Unmerged, Uncorrected> {}
impl CanCorrectPairUnchecked for CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected> {}
impl CanCorrectPairUnchecked for CanTuple<HasOverlap, Validated, Unmerged, Uncorrected> {}
impl CanCorrectMerged for MergedRead {
    fn into_corrected_merged(self) -> Result<CorrectedMergedRead> {
        self.into_uncorrected().correct()
    }
}

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

pub(crate) trait CorrectPairOp {
    type Out;
    fn correct_pair(self) -> Result<Self::Out>;
}

pub(crate) trait CorrectPairUncheckedOp {
    type Out;
    fn correct_pair_unchecked(self) -> Result<Self::Out>;
}

pub(crate) trait CorrectMergedOp {
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

        let overlap_outcome = match read_pair.overlap(&assembler.config().overlap)? {
            Some(overlap) => OverlapOutcome::Found(OverlapSnapshot::from_overlap(&overlap)),
            None => OverlapOutcome::Missing,
        };

        Ok(PairContext {
            assembler,
            input,
            read_pair,
            overlap_outcome,
            _marker: PhantomData,
        })
    }
}

impl<'asm, 'pair, R, O, V, M, C> ValidateOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanValidate,
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap + HasReadPair,
{
    type Out = ValidatedContext<'asm, 'pair, R>;

    fn validate(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(ctx.read_pair_ref());
            overlap.validate(ctx.read_pair_ref(), &ctx.assembler_ref().config().validator)?;
            let (assembler, input, read_pair) = ctx.into_parts();

            Ok(PairContext {
                assembler,
                input,
                read_pair,
                overlap_outcome: OverlapOutcome::Found(snapshot),
                _marker: PhantomData,
            })
        })?
        .on_missing(|ctx| {
            let (assembler, input, read_pair) = ctx.into_parts();
            Ok(PairContext {
                assembler,
                input,
                read_pair,
                overlap_outcome: OverlapOutcome::Missing,
                _marker: PhantomData,
            })
        })
    }
}

impl<'asm, 'pair, R, O, V, M, C> ValidatePredicateOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap + HasReadPair,
{
    fn is_valid(&self) -> Result<bool> {
        match self.overlap_outcome() {
            OverlapOutcome::Found(snapshot) => {
                let overlap = snapshot.materialize_overlap(self.read_pair_ref());
                Ok(overlap
                    .validate(
                        self.read_pair_ref(),
                        &self.assembler_ref().config().validator,
                    )
                    .is_ok())
            },
            OverlapOutcome::Missing | OverlapOutcome::Unknown => Ok(false),
        }
    }
}

impl<'asm, 'pair, R> MergeOp
    for PairContext<'asm, 'pair, R, HasOverlap, Unvalidated, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected>: CanMerge,
    PairContext<'asm, 'pair, R, HasOverlap, Unvalidated, Unmerged, Uncorrected>:
        HasPairOverlap + HasReadPair,
{
    type Out = MergedRead;

    fn merge(self) -> Result<Self::Out> {
        self.on_found(|ctx, _snapshot| merge_from(&ctx))?
            .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<'asm, 'pair, R> MergeOp
    for PairContext<'asm, 'pair, R, HasOverlap, Validated, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, Validated, Unmerged, Uncorrected>: CanMerge,
    PairContext<'asm, 'pair, R, HasOverlap, Validated, Unmerged, Uncorrected>:
        HasPairOverlap + HasReadPair,
{
    type Out = MergedRead;

    fn merge(self) -> Result<Self::Out> {
        self.on_found(|ctx, _snapshot| merge_from(&ctx))?
            .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<'asm, 'pair, R, O, V, M, C> CorrectPairOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanCorrectPair,
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap + HasReadPair,
{
    type Out = CorrectedReadPair;

    fn correct_pair(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(ctx.read_pair_ref());
            let validated =
                overlap.validate(ctx.read_pair_ref(), &ctx.assembler_ref().config().validator)?;
            ctx.read_pair_ref()
                .correct_from_overlap(validated.overlap())
        })?
        .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<'asm, 'pair, R, O, V, M, C> CorrectPairUncheckedOp for PairContext<'asm, 'pair, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanCorrectPairUnchecked,
    PairContext<'asm, 'pair, R, O, V, M, C>: HasPairOverlap + HasReadPair,
{
    type Out = CorrectedReadPair;

    fn correct_pair_unchecked(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(ctx.read_pair_ref());
            ctx.read_pair_ref().correct_from_overlap(&overlap)
        })?
        .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<T> CorrectMergedOp for T
where
    T: CanCorrectMerged,
{
    type Out = CorrectedMergedRead;

    fn correct(self) -> Result<Self::Out> {
        self.into_corrected_merged()
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
    pub fn process(self) -> Result<CorrectedMergedRead> {
        let merged = MergeOp::merge(ValidateOp::validate(OverlapOp::overlap(self)?)?)?;
        CorrectMergedOp::correct(merged)
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

    /// Merge this overlap without validation checks.
    ///
    /// # Errors
    ///
    /// Returns an error if merge fails.
    pub fn merge_unchecked(self) -> Result<MergedRead> {
        MergeOp::merge(self)
    }

    /// Correct both mates directly from overlap evidence without validation.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap-derived correction fails.
    pub fn correct_pair_unchecked(self) -> Result<CorrectedReadPair> {
        CorrectPairUncheckedOp::correct_pair_unchecked(self)
    }
}

impl<R> ValidatedContext<'_, '_, R>
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
    pub fn merge(self) -> Result<MergedRead> {
        MergeOp::merge(self)
    }

    /// Correct both mates directly from overlap evidence after validation.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap-derived correction fails.
    pub fn correct_pair(self) -> Result<CorrectedReadPair> {
        CorrectPairOp::correct_pair(self)
    }
}
