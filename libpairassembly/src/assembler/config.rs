//! Top-level assembler configuration and orchestration types.

use std::marker::PhantomData;

use crate::{
    BaseCallValidator, OverlapParams, Result,
    correct::{CorrectedMergedRead, CorrectionParams},
};

use super::{
    PairContext, PairInput, PairReady, ProcessIter, SeqRecordView, context::OverlapOutcome,
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
    pub fn process_pair<R>(&self, pair: &PairInput<R>) -> Result<CorrectedMergedRead>
    where
        R: SeqRecordView,
    {
        Ok(self
            .on_pair(pair)?
            .overlap()?
            .validate()?
            .merge()?
            .correct()?
            .into_corrected_merged_read())
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
            validation_metrics: None,
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
