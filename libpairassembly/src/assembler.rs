//! Assembler-centered public API with two usage layers.
//!
//! - Layer A: per-pair fluent DAG transitions (`on_pair(...).overlap()...`).
//! - Layer B: collection orchestration (`process_iter`, `process_iter_with`).
//!
//! The default convenience path (`process_pair` / `process_iter`) is the checked
//! path (`overlap -> validate -> merge -> correct`). Unchecked expert paths are
//! available from [`OverlapContext`] (`merge_unchecked`,
//! `correct_pair_unchecked`) and are intentionally explicit.
//!
//! Transition channels are tracked across four dimensions:
//! - `O`: overlap discovered (`NoOverlap`/`HasOverlap`)
//! - `V`: overlap validated (`Unvalidated`/`Validated`)
//! - `M`: merged state (`Unmerged`/`Merged`)
//! - `C`: correction state (`Uncorrected`/`PairCorrected`/`MergedCorrected`)
//!
//! ```text
//! PairReady
//!   O=NoOverlap, V=Unvalidated, M=Unmerged, C=Uncorrected
//!       |
//!       | overlap()
//!       v
//! OverlapContext
//!   O=HasOverlap, V=Unvalidated, M=Unmerged, C=Uncorrected
//!     |                    |\
//!     | validate()         | merge_unchecked() -> correct()
//!     |                    |  \
//!     v                    |   > CorrectedMergedRead (unchecked branch)
//! ValidatedContext         |
//!   O=HasOverlap,          | correct_pair_unchecked()
//!   V=Validated,           v
//!   M=Unmerged,            CorrectedReadPair (unchecked branch)
//!   C=Uncorrected
//!     |            |
//!     | merge()    | correct_pair()
//!     v            v
//! correct()        CorrectedReadPair
//!     |
//!     v
//! CorrectedMergedRead (checked branch)
//! ```

use std::{fmt::Display, marker::PhantomData};

use crate::{
    BaseCallValidator, OverlapParams, PairOverlap, ReadPair, Result, SequenceRead, TiePolicy,
    correct::{CorrectedMergedRead, CorrectedReadPair, CorrectionParams},
    errors::{ConversionError, OverlapError},
    merge::UncorrectedMergedRead,
};
use state::{
    HasOverlap, Merged, MergedCorrected, NoOverlap, PairCorrected, Uncorrected, Unmerged,
    Unvalidated, Validated,
};

/// Placeholder merge-stage configuration for the top-level `Assembler` API.
///
/// This currently carries no options and exists to preserve API shape while merge
/// behavior is further decomposed.
#[derive(Debug, Clone, Copy, Default)]
pub struct MergeParams;

/// Execution strategy for collection processing entrypoints.
///
/// - [`ExecutionPolicy::Record`] processes each pair independently.
/// - [`ExecutionPolicy::Batch`] reserves an explicit policy surface for future
///   data-oriented batch execution.
#[derive(Debug, Clone, Copy, Default)]
pub enum ExecutionPolicy {
    #[default]
    Record,
    Batch(BatchPolicy),
}

impl ExecutionPolicy {
    /// Select per-pair record execution.
    #[must_use]
    pub fn record() -> Self {
        Self::Record
    }

    /// Select batch execution with default batch policy values.
    #[must_use]
    pub fn batch() -> Self {
        Self::Batch(BatchPolicy::default())
    }
}

/// Configuration knobs for batch execution.
///
/// These fields are currently part of API scaffolding and may gain stricter
/// semantics as batch backend implementation matures.
#[derive(Debug, Clone, Copy)]
pub struct BatchPolicy {
    pub chunk_pairs: usize,
    pub max_bytes: usize,
    pub precompute_revcomp: bool,
    pub threads: Option<usize>,
}

impl Default for BatchPolicy {
    fn default() -> Self {
        Self {
            chunk_pairs: 1024,
            max_bytes: 64 * 1024 * 1024,
            precompute_revcomp: true,
            threads: None,
        }
    }
}

/// Top-level assembler configuration.
///
/// This bundles stage-specific settings and execution strategy in one place so
/// callers can configure and reuse an `Assembler` across many pairs.
#[derive(Debug, Clone)]
pub struct AssemblerConfig {
    pub overlap: OverlapParams,
    pub validator: BaseCallValidator,
    pub merge: MergeParams,
    pub correction: CorrectionParams,
    pub execution: ExecutionPolicy,
}

impl Default for AssemblerConfig {
    fn default() -> Self {
        Self {
            overlap: OverlapParams::default(),
            validator: BaseCallValidator::default(),
            merge: MergeParams,
            correction: CorrectionParams::default(),
            execution: ExecutionPolicy::default(),
        }
    }
}

/// Top-level API object for pair assembly orchestration.
///
/// `Assembler` currently delegates to existing overlap/validate/merge/correct
/// internals while exposing a stable surface for the in-progress API migration.
#[derive(Debug, Clone)]
pub struct Assembler {
    config: AssemblerConfig,
}

impl Assembler {
    /// Start building a configured [`Assembler`].
    #[must_use]
    pub fn builder() -> AssemblerBuilder {
        AssemblerBuilder::default()
    }

    /// Borrow the active configuration.
    #[must_use]
    pub fn config(&self) -> &AssemblerConfig {
        &self.config
    }

    /// Process a single paired input record to a corrected merged read.
    ///
    /// This convenience method runs the canonical checked path:
    /// `overlap -> validate -> merge -> correct`.
    ///
    /// # Errors
    ///
    /// Returns an error if pairing, overlap discovery, validation, merging, or
    /// correction fail for this input pair.
    pub fn process_pair<R>(&self, pair: PairInput<R>) -> Result<CorrectedMergedRead>
    where
        R: SeqRecordView,
    {
        self.on_pair(&pair)?
            .overlap()?
            .validate()?
            .merge()?
            .correct()
    }

    /// Process an iterator of paired records with this assembler configuration.
    ///
    /// Each output item corresponds to one input pair and is returned as a
    /// `Result` so callers can decide whether to fail-fast or handle per-pair
    /// errors inline.
    pub fn process_iter<'asm, I, R>(&'asm self, pairs: I) -> ProcessIter<'asm, I::IntoIter>
    where
        I: IntoIterator<Item = PairInput<R>> + 'asm,
        R: SeqRecordView + 'asm,
    {
        ProcessIter {
            assembler: self,
            iter: pairs.into_iter(),
            execution: self.config.execution,
        }
    }

    /// Process an iterator of paired records with a custom per-pair pipeline closure.
    ///
    /// This advanced entrypoint enables callers to choose explicit branch ordering
    /// for each pair while preserving the same per-item `Result` behavior as
    /// [`Assembler::process_iter`]. Use this when you intentionally need unchecked
    /// expert branches or mixed terminal output types.
    pub fn process_iter_with<'asm, I, R, O, F>(
        &'asm self,
        pairs: I,
        mut f: F,
    ) -> impl Iterator<Item = Result<O>> + 'asm
    where
        I: IntoIterator<Item = PairInput<R>> + 'asm,
        R: SeqRecordView + 'asm,
        F: for<'pair> FnMut(PairReady<'asm, 'pair, R>) -> Result<O> + 'asm,
    {
        pairs.into_iter().map(move |pair| {
            let ready = self.on_pair(&pair)?;
            f(ready)
        })
    }

    /// Begin explicit per-pair processing flow.
    ///
    /// This entrypoint exists to preserve fluent per-pair APIs while collection
    /// APIs are layered on top.
    ///
    /// # Errors
    ///
    /// Returns an error only if pair initialization fails.
    pub fn on_pair<'pair, R>(&self, pair: &'pair PairInput<R>) -> Result<PairReady<'_, 'pair, R>>
    where
        R: SeqRecordView,
    {
        let read_pair = pair.try_into_read_pair()?;
        Ok(PairContext {
            assembler: self,
            input: pair,
            read_pair,
            overlap_outcome: OverlapOutcome::Unknown,
            _marker: PhantomData,
        })
    }
}

/// Builder for [`Assembler`].
#[derive(Debug, Clone, Default)]
pub struct AssemblerBuilder {
    config: AssemblerConfig,
}

impl AssemblerBuilder {
    /// Set overlap detection parameters.
    #[must_use]
    pub fn overlap(mut self, overlap: OverlapParams) -> Self {
        self.config.overlap = overlap;
        self
    }

    /// Set overlap validation parameters.
    #[must_use]
    pub fn validate(mut self, validator: BaseCallValidator) -> Self {
        self.config.validator = validator;
        self
    }

    /// Set merge-stage parameters.
    #[must_use]
    pub fn merge(mut self, merge: MergeParams) -> Self {
        self.config.merge = merge;
        self
    }

    /// Set quality-correction parameters.
    #[must_use]
    pub fn correct(mut self, correction: CorrectionParams) -> Self {
        self.config.correction = correction;
        self
    }

    /// Set execution policy.
    #[must_use]
    pub fn execution(mut self, execution: ExecutionPolicy) -> Self {
        self.config.execution = execution;
        self
    }

    /// Finalize and return the configured assembler.
    ///
    /// # Errors
    ///
    /// This currently does not fail in practice, but returns `Result` to keep
    /// construction compatible with future validation during build.
    pub fn build(self) -> Result<Assembler> {
        Ok(Assembler {
            config: self.config,
        })
    }
}

mod state {
    #[derive(Debug, Clone, Copy)]
    pub struct NoOverlap;
    #[derive(Debug, Clone, Copy)]
    pub struct HasOverlap;

    #[derive(Debug, Clone, Copy)]
    pub struct Unvalidated;
    #[derive(Debug, Clone, Copy)]
    pub struct Validated;

    #[derive(Debug, Clone, Copy)]
    pub struct Unmerged;
    #[derive(Debug, Clone, Copy)]
    pub struct Merged;

    #[derive(Debug, Clone, Copy)]
    pub struct Uncorrected;
    #[derive(Debug, Clone, Copy)]
    pub struct PairCorrected;
    #[derive(Debug, Clone, Copy)]
    pub struct MergedCorrected;
}

/// Internal typestate carrier for per-pair Assembler DAG transitions.
#[derive(Debug, Clone)]
pub struct PairContext<'asm, 'pair, R, O, V, M, C> {
    assembler: &'asm Assembler,
    input: &'pair PairInput<R>,
    read_pair: ReadPair<'pair>,
    overlap_outcome: OverlapOutcome,
    _marker: PhantomData<(O, V, M, C)>,
}

#[derive(Debug, Clone, Copy)]
enum OverlapOutcome {
    Unknown,
    Missing,
    Found(OverlapSnapshot),
}

#[derive(Debug)]
enum OverlapBranch<C, T> {
    Value(T),
    Context(C),
}

impl<C, T> OverlapBranch<C, T> {
    fn on_missing(self, f: impl FnOnce(C) -> Result<T>) -> Result<T> {
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
impl CanCorrectMerged for UncorrectedMergedRead {
    fn into_corrected_merged(self) -> Result<CorrectedMergedRead> {
        UncorrectedMergedRead::correct(self)
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

/// Iterator adaptor for processing paired records through an [`Assembler`].
#[derive(Debug)]
pub struct ProcessIter<'asm, I> {
    assembler: &'asm Assembler,
    iter: I,
    execution: ExecutionPolicy,
}

impl<I, R> Iterator for ProcessIter<'_, I>
where
    I: Iterator<Item = PairInput<R>>,
    R: SeqRecordView,
{
    type Item = Result<CorrectedMergedRead>;

    fn next(&mut self) -> Option<Self::Item> {
        let pair = self.iter.next()?;
        let result = match self.execution {
            ExecutionPolicy::Record => self.assembler.process_pair(pair),
            ExecutionPolicy::Batch(_policy) => self.assembler.process_pair(pair),
        };
        Some(result)
    }
}

impl<R, O, V, M, C> PairContext<'_, '_, R, O, V, M, C> {
    #[inline]
    fn on_found<T>(
        self,
        f: impl FnOnce(Self, OverlapSnapshot) -> Result<T>,
    ) -> Result<OverlapBranch<Self, T>> {
        match self.overlap_outcome {
            OverlapOutcome::Found(snapshot) => Ok(OverlapBranch::Value(f(self, snapshot)?)),
            OverlapOutcome::Missing | OverlapOutcome::Unknown => Ok(OverlapBranch::Context(self)),
        }
    }

    #[inline]
    fn on_missing<T>(self, f: impl FnOnce(Self) -> Result<T>) -> Result<OverlapBranch<Self, T>> {
        match self.overlap_outcome {
            OverlapOutcome::Missing => Ok(OverlapBranch::Value(f(self)?)),
            OverlapOutcome::Found(_) | OverlapOutcome::Unknown => Ok(OverlapBranch::Context(self)),
        }
    }
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

        let overlap_outcome = match read_pair.overlap(&assembler.config.overlap)? {
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
{
    type Out = ValidatedContext<'asm, 'pair, R>;

    fn validate(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(&ctx.read_pair);
            overlap.validate(&ctx.read_pair, &ctx.assembler.config.validator)?;

            Ok(PairContext {
                assembler: ctx.assembler,
                input: ctx.input,
                read_pair: ctx.read_pair,
                overlap_outcome: OverlapOutcome::Found(snapshot),
                _marker: PhantomData,
            })
        })?
        .on_missing(|ctx| {
            Ok(PairContext {
                assembler: ctx.assembler,
                input: ctx.input,
                read_pair: ctx.read_pair,
                overlap_outcome: OverlapOutcome::Missing,
                _marker: PhantomData,
            })
        })
    }
}

impl<R> MergeOp for PairContext<'_, '_, R, HasOverlap, Unvalidated, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, Unvalidated, Unmerged, Uncorrected>: CanMerge,
{
    type Out = UncorrectedMergedRead;

    fn merge(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(&ctx.read_pair);
            let validated = crate::ValidatedOverlap {
                mates: &ctx.read_pair,
                overlap,
            };
            validated.merge()
        })?
        .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<R> MergeOp for PairContext<'_, '_, R, HasOverlap, Validated, Unmerged, Uncorrected>
where
    R: SeqRecordView,
    CanTuple<HasOverlap, Validated, Unmerged, Uncorrected>: CanMerge,
{
    type Out = UncorrectedMergedRead;

    fn merge(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(&ctx.read_pair);
            overlap
                .validate(&ctx.read_pair, &ctx.assembler.config.validator)?
                .merge()
        })?
        .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<R, O, V, M, C> CorrectPairOp for PairContext<'_, '_, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanCorrectPair,
{
    type Out = CorrectedReadPair;

    fn correct_pair(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(&ctx.read_pair);
            let validated = overlap.validate(&ctx.read_pair, &ctx.assembler.config.validator)?;
            ctx.read_pair.correct_from_overlap(&validated.overlap)
        })?
        .on_missing(|_| Err(OverlapError::NoOverlapFound.into()))
    }
}

impl<R, O, V, M, C> CorrectPairUncheckedOp for PairContext<'_, '_, R, O, V, M, C>
where
    R: SeqRecordView,
    CanTuple<O, V, M, C>: CanCorrectPairUnchecked,
{
    type Out = CorrectedReadPair;

    fn correct_pair_unchecked(self) -> Result<Self::Out> {
        self.on_found(|ctx, snapshot| {
            let overlap = snapshot.materialize_overlap(&ctx.read_pair);
            ctx.read_pair.correct_from_overlap(&overlap)
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
    pub fn merge_unchecked(self) -> Result<UncorrectedMergedRead> {
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
    /// Merge this pair using the checked path.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation or merge fails.
    pub fn merge(self) -> Result<UncorrectedMergedRead> {
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

#[derive(Debug, Clone, Copy)]
struct OverlapSnapshot {
    overlap_len: usize,
    r1_start_offset: usize,
    r1_end_offset: usize,
    r2_start_offset: usize,
    r2_end_offset: usize,
}

impl OverlapSnapshot {
    fn from_overlap(overlap: &PairOverlap<'_>) -> Self {
        Self {
            overlap_len: overlap.len(),
            r1_start_offset: overlap.forward_start_offset(),
            r1_end_offset: overlap.forward_end_offset(),
            r2_start_offset: overlap.reverse_start_offset(),
            r2_end_offset: overlap.reverse_end_offset(),
        }
    }

    fn materialize_overlap<'pair>(self, pair: &'pair ReadPair<'pair>) -> PairOverlap<'pair> {
        let r1_seq = pair.fwd_mate.sequence().as_bytes();
        let r1_qual = pair.fwd_mate.quality_scores().as_bytes();
        let r2_seq_rc = pair.rev_mate.reverse_complement();
        let mut r2_qual_rc = pair.rev_mate.quality_scores().as_bytes().to_vec();
        r2_qual_rc.reverse();

        PairOverlap::from_components(
            self.overlap_len,
            self.r1_start_offset,
            self.r1_end_offset,
            self.r2_start_offset,
            self.r2_end_offset,
            &r1_seq[self.r1_start_offset..=self.r1_end_offset],
            &r1_qual[self.r1_start_offset..=self.r1_end_offset],
            r2_seq_rc[self.r2_start_offset..=self.r2_end_offset].to_vec(),
            r2_qual_rc[self.r2_start_offset..=self.r2_end_offset].to_vec(),
        )
    }
}

/// Boundary trait for pair records accepted by the assembler API.
///
/// Implement this for external record types to use `Assembler` directly.
pub trait SeqRecordView {
    fn id(&self) -> &str;
    fn seq(&self) -> &str;
    fn qual(&self) -> &str;
}

/// Boundary trait for constructing user-space record types from corrected output parts.
pub trait FromRecordParts: Sized {
    type Error;

    /// Construct a record instance from owned identifier, sequence, and quality parts.
    ///
    /// # Errors
    ///
    /// Returns an error if the target record type cannot be constructed from the
    /// provided parts.
    fn try_from_parts(
        id: String,
        seq: Vec<u8>,
        qual: Vec<u8>,
    ) -> std::result::Result<Self, Self::Error>;
}

/// Trait for extracting owned merged-record parts from terminal pipeline outputs.
pub trait IntoOwnedRecordParts {
    fn into_owned_record_parts(self) -> (String, Vec<u8>, Vec<u8>);
}

/// Trait for extracting owned paired-record parts from terminal pipeline outputs.
pub trait IntoOwnedPairRecordParts {
    fn into_owned_pair_record_parts(self) -> (String, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);
}

/// Blanket conversion helper for merged outputs.
///
/// By convention, `into_*` methods in this crate are reserved for meaningful
/// conversion across representation boundaries. Identity-shaped `into_*`
/// endpoints are intentionally not provided.
pub trait IntoRecordConversion: IntoOwnedRecordParts {
    /// Convert merged terminal output into a user record type.
    ///
    /// # Errors
    ///
    /// Returns an error if the target record type cannot be constructed from
    /// the owned parts extracted from `self`.
    fn into_record<T>(self) -> Result<T>
    where
        Self: Sized,
        T: FromRecordParts,
        T::Error: Display,
    {
        let (id, seq, qual) = self.into_owned_record_parts();
        T::try_from_parts(id, seq, qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()).into())
    }
}

impl<T> IntoRecordConversion for T where T: IntoOwnedRecordParts {}

/// Blanket conversion helper for paired outputs.
///
/// By convention, `into_*` methods in this crate are reserved for meaningful
/// conversion across representation boundaries. Identity-shaped `into_*`
/// endpoints are intentionally not provided.
pub trait IntoRecordsConversion: IntoOwnedPairRecordParts {
    /// Convert paired terminal output into two user record values.
    ///
    /// # Errors
    ///
    /// Returns an error if either target record value cannot be constructed
    /// from the owned parts extracted from `self`.
    fn into_records<T>(self) -> Result<(T, T)>
    where
        Self: Sized,
        T: FromRecordParts,
        T::Error: Display,
    {
        let (id, fwd_seq, fwd_qual, rev_seq, rev_qual) = self.into_owned_pair_record_parts();
        let left = T::try_from_parts(id.clone(), fwd_seq, fwd_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        let right = T::try_from_parts(id, rev_seq, rev_qual)
            .map_err(|err| ConversionError::RecordConstruction(err.to_string()))?;
        Ok((left, right))
    }
}

impl<T> IntoRecordsConversion for T where T: IntoOwnedPairRecordParts {}

impl SeqRecordView for SequenceRead<'_> {
    fn id(&self) -> &str {
        self.id()
    }

    fn seq(&self) -> &str {
        self.sequence()
    }

    fn qual(&self) -> &str {
        self.quality_scores()
    }
}

/// Pair wrapper accepted by assembler entrypoints.
#[derive(Debug)]
pub struct PairInput<R> {
    pub r1: R,
    pub r2: R,
}

impl<R> PairInput<R> {
    /// Construct a paired input wrapper.
    #[must_use]
    pub fn new(r1: R, r2: R) -> Self {
        Self { r1, r2 }
    }

    /// Convert a generic pair input into the canonical internal [`ReadPair`] form.
    ///
    /// # Errors
    ///
    /// Returns an error if either record has invalid sequence/quality structure or
    /// if IDs do not correspond to a valid read pair.
    pub fn try_into_read_pair(&self) -> Result<ReadPair<'_>>
    where
        R: SeqRecordView,
    {
        let read1 = SequenceRead::try_new(self.r1.id(), self.r1.seq(), self.r1.qual())?;
        let read2 = SequenceRead::try_new(self.r2.id(), self.r2.seq(), self.r2.qual())?;
        ReadPair::from(read1, read2)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;
    use crate::{Error, errors::OverlapError, test_fixtures::TupleRecord};

    fn rec(id: &str, seq: &str, qual: &str) -> TupleRecord {
        TupleRecord::from_strs(id, seq, qual)
    }

    fn demo_pair(id: &str) -> PairInput<TupleRecord> {
        PairInput::new(rec(id, "TTTACGTA", "IIIIIIII"), rec(id, "TACGT", "IIIII"))
    }

    #[test]
    fn test_builder_with_defaults() {
        let asm = Assembler::builder()
            .build()
            .expect("default assembler builder should produce a valid configuration");
        assert!(matches!(asm.config().execution, ExecutionPolicy::Record));
    }

    #[test]
    fn test_process_pair_with_tuple_record_fixture() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read1", "ACGTACGT", "IIIIIIII"),
            rec("read1", "ACGTACGT", "IIIIIIII"),
        );

        let result = asm.process_pair(pair);
        assert!(matches!(
            result,
            Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
        ));
    }

    #[test]
    fn test_process_iter_yields_results() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pairs = vec![demo_pair("read1"), demo_pair("read2")];

        let results = asm.process_iter(pairs).collect::<Vec<_>>();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_on_pair_process_delegates() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read1", "ACGTACGT", "IIIIIIII"),
            rec("read1", "ACGTACGT", "IIIIIIII"),
        );

        let delegated = asm
            .on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .process();
        assert!(matches!(
            delegated,
            Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
        ));
    }

    #[test]
    fn test_context_checked_and_unchecked_paths_exist() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair1 = demo_pair("read1");
        let pair2 = demo_pair("read2");

        let checked = asm
            .on_pair(&pair1)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors")
            .validate();
        assert!(checked.is_ok());

        let unchecked = asm
            .on_pair(&pair2)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors")
            .merge_unchecked();
        assert!(unchecked.is_ok());
    }

    #[test]
    fn test_process_pair_equals_process_iter_singleton_success() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");

        let single = asm
            .process_pair(demo_pair("read-single"))
            .expect("singleton process_pair should succeed for demo pair");
        let iter = asm
            .process_iter(vec![demo_pair("read-single")])
            .next()
            .expect("iterator should yield one singleton result")
            .expect("singleton process_iter result should succeed for demo pair");

        assert_eq!(single.id(), iter.id());
        assert_eq!(single.sequence_bytes(), iter.sequence_bytes());
        assert_eq!(single.quality_bytes(), iter.quality_bytes());
    }

    #[test]
    fn test_process_pair_equals_process_iter_singleton_error() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read-tie", "ACGTACGT", "IIIIIIII"),
            rec("read-tie", "ACGTACGT", "IIIIIIII"),
        );

        let single = asm.process_pair(pair).unwrap_err();
        assert!(matches!(
            single,
            Error::OverlapError(OverlapError::OverlapTie { .. })
        ));

        let iter = asm
            .process_iter(vec![PairInput::new(
                rec("read-tie", "ACGTACGT", "IIIIIIII"),
                rec("read-tie", "ACGTACGT", "IIIIIIII"),
            )])
            .next()
            .expect("iterator should yield one singleton error result")
            .unwrap_err();
        assert!(matches!(
            iter,
            Error::OverlapError(OverlapError::OverlapTie { .. })
        ));
    }

    #[test]
    fn test_process_iter_batch_policy_matches_record_policy() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm_record = Assembler::builder()
            .overlap(overlap)
            .execution(ExecutionPolicy::record())
            .build()
            .expect("record execution policy should build successfully");
        let asm_batch = Assembler::builder()
            .overlap(overlap)
            .execution(ExecutionPolicy::batch())
            .build()
            .expect("batch execution policy should build successfully");

        let record = asm_record
            .process_iter(vec![demo_pair("read-policy")])
            .next()
            .expect("record-policy iterator should yield a singleton result")
            .expect("record-policy singleton result should succeed");
        let batch = asm_batch
            .process_iter(vec![demo_pair("read-policy")])
            .next()
            .expect("batch-policy iterator should yield a singleton result")
            .expect("batch-policy singleton result should succeed");

        assert_eq!(record.id(), batch.id());
        assert_eq!(record.sequence_bytes(), batch.sequence_bytes());
        assert_eq!(record.quality_bytes(), batch.quality_bytes());
    }

    #[test]
    fn test_overlap_context_clone_branches_without_recomputing_overlap_selection() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = demo_pair("read-clone");

        let ctx = asm
            .on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors");
        let checked = ctx
            .clone()
            .validate()
            .expect("validation should succeed for overlap-clone fixture")
            .merge()
            .expect("checked merge should succeed for overlap-clone fixture")
            .correct()
            .expect("checked correction should succeed for overlap-clone fixture");
        let unchecked = ctx
            .merge_unchecked()
            .expect("unchecked merge should succeed for overlap-clone fixture")
            .correct()
            .expect("unchecked correction should succeed for overlap-clone fixture");

        assert_eq!(checked.id(), unchecked.id());
        assert_eq!(checked.sequence_bytes(), unchecked.sequence_bytes());
        assert_eq!(checked.quality_bytes(), unchecked.quality_bytes());
    }

    #[test]
    fn test_correct_pair_checked_and_unchecked_paths_match() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = demo_pair("read-correct");

        let ctx = asm
            .on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors");
        let checked = ctx
            .clone()
            .validate()
            .expect("validation should succeed for checked-vs-unchecked fixture")
            .correct_pair()
            .expect("checked correction should succeed for checked-vs-unchecked fixture");
        let unchecked = ctx
            .correct_pair_unchecked()
            .expect("unchecked correction should succeed for checked-vs-unchecked fixture");

        assert_eq!(checked.id(), unchecked.id());
        assert_eq!(checked.fwd_sequence_bytes(), unchecked.fwd_sequence_bytes());
        assert_eq!(checked.fwd_quality_bytes(), unchecked.fwd_quality_bytes());
        assert_eq!(checked.rev_sequence_bytes(), unchecked.rev_sequence_bytes());
        assert_eq!(checked.rev_quality_bytes(), unchecked.rev_quality_bytes());
    }

    #[test]
    fn test_correct_pair_checked_path_fails_for_low_confidence_overlap() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let validator = BaseCallValidator::new().with_min_entropy(44);
        let asm = Assembler::builder()
            .overlap(overlap)
            .validate(validator)
            .build()
            .expect("assembler builder should accept explicit overlap/validation settings");
        let pair = PairInput::new(
            rec("read-low-confidence", "ACGTACGT", "IIIIIIII"),
            rec("read-low-confidence", "TCGTACGT", "IIIIIIII"),
        );

        let ctx = asm
            .on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors");
        assert!(ctx.clone().correct_pair_unchecked().is_ok());
        assert!(
            ctx.validate()
                .and_then(ValidatedContext::correct_pair)
                .is_err()
        );
    }

    #[test]
    fn test_no_overlap_outcome_flows_through_context_and_fails_on_consumers() {
        let overlap = OverlapParams::default()
            .with_min_overlap(4)
            .with_min_comparisons(4);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read-no-overlap", "AAAAAAAA", "IIIIIIII"),
            rec("read-no-overlap", "CCCCCCCC", "IIIIIIII"),
        );

        let overlapped = asm
            .on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap()
            .expect("overlap stage should run without scanner/conversion errors");
        let validated = overlapped
            .clone()
            .validate()
            .expect("validation should succeed even when no-overlap is carried forward");

        assert!(matches!(
            overlapped.clone().merge_unchecked(),
            Err(Error::OverlapError(OverlapError::NoOverlapFound))
        ));
        assert!(matches!(
            validated.clone().merge(),
            Err(Error::OverlapError(OverlapError::NoOverlapFound))
        ));
        assert!(matches!(
            overlapped.correct_pair_unchecked(),
            Err(Error::OverlapError(OverlapError::NoOverlapFound))
        ));
        assert!(matches!(
            validated.correct_pair(),
            Err(Error::OverlapError(OverlapError::NoOverlapFound))
        ));
    }

    #[test]
    fn test_process_pair_reports_no_overlap_outcome_at_merge_stage() {
        let overlap = OverlapParams::default()
            .with_min_overlap(4)
            .with_min_comparisons(4);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read-no-overlap-process", "AAAAAAAA", "IIIIIIII"),
            rec("read-no-overlap-process", "CCCCCCCC", "IIIIIIII"),
        );

        assert!(matches!(
            asm.process_pair(pair),
            Err(Error::OverlapError(OverlapError::NoOverlapFound))
        ));
    }

    #[test]
    fn test_process_iter_singleton_no_overlap_matches_process_pair_error() {
        let overlap = OverlapParams::default()
            .with_min_overlap(4)
            .with_min_comparisons(4);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
            rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
        );

        let single = asm.process_pair(pair).unwrap_err();
        let iter = asm
            .process_iter(vec![PairInput::new(
                rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
                rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
            )])
            .next()
            .expect("iterator should yield one singleton error result")
            .unwrap_err();

        assert!(matches!(
            single,
            Error::OverlapError(OverlapError::NoOverlapFound)
        ));
        assert!(matches!(
            iter,
            Error::OverlapError(OverlapError::NoOverlapFound)
        ));
    }

    #[test]
    fn test_overlap_tie_still_errors_at_overlap_stage() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pair = PairInput::new(
            rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
            rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
        );

        assert!(matches!(
            asm.on_pair(&pair)
                .expect("on_pair should convert tuple records into read-pair context")
                .overlap(),
            Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
        ));
    }

    #[test]
    fn test_process_iter_with_custom_checked_merge_pipeline() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pairs = vec![demo_pair("read-custom-1"), demo_pair("read-custom-2")];

        let results = asm
            .process_iter_with(pairs, |ready| {
                ready.overlap()?.validate()?.merge()?.correct()
            })
            .collect::<Vec<_>>();

        assert_eq!(results.len(), 2);
        assert!(results.iter().all(Result::is_ok));
    }

    #[test]
    fn test_process_iter_with_custom_unmerged_pipeline() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder()
            .overlap(overlap)
            .build()
            .expect("assembler builder should accept explicit overlap settings");
        let pairs = vec![demo_pair("read-custom-unmerged")];

        let result = asm
            .process_iter_with(pairs, |ready| ready.overlap()?.correct_pair_unchecked())
            .next()
            .expect("iterator should yield one singleton custom-pipeline result")
            .expect("custom unchecked pipeline should succeed for demo pair");

        assert_eq!(result.id(), "read-custom-unmerged");
    }
}
