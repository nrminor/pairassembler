//! Top-level assembler configuration and orchestration types.

use std::marker::PhantomData;

use crate::{
    OverlapParams, OverlapValidator, Result, correct::CorrectionParams, merge::MergeParams,
    read::OwnedSequenceRead,
};

use super::{PairInput, PairReady, ProcessIter, SeqRecordView, context::PairContext};

/// Top-level assembler configuration.
///
/// This bundles stage-specific settings in one place so callers can configure and reuse an
/// [`Assembler`] across many pairs.
#[derive(Debug, Clone, Default)]
pub struct AssemblerConfig {
    pub overlap: OverlapParams,
    pub validator: OverlapValidator,
    pub merge: MergeParams,
    pub correction: CorrectionParams,
}

/// Top-level API object for pair assembly orchestration.
///
/// `Assembler` is cheap to clone and intended to be configured once, then reused across many input
/// pairs.
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

    #[inline]
    pub(crate) fn overlap_params(&self) -> &OverlapParams {
        &self.config.overlap
    }

    #[inline]
    pub(crate) fn validator(&self) -> &OverlapValidator {
        &self.config.validator
    }

    #[inline]
    pub(crate) fn correction_params(&self) -> CorrectionParams {
        self.config.correction
    }

    #[inline]
    pub(crate) fn merge_params(&self) -> MergeParams {
        self.config.merge
    }

    /// Process a single paired input record to a corrected merged read when an overlap is found.
    ///
    /// This convenience method runs the canonical checked path:
    /// `find_overlap -> validate -> merge -> correct`.
    ///
    /// ```rust
    /// use libpairassembly::prelude::*;
    ///
    /// # fn main() -> libpairassembly::Result<()> {
    /// let pair = PairInput::new(
    ///     SequenceRead::try_new(
    ///         "read-1",
    ///         "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
    ///         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ///     )?,
    ///     SequenceRead::try_new(
    ///         "read-1",
    ///         "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
    ///         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ///     )?,
    /// );
    ///
    /// let merged = Assembler::builder()
    ///     .build()?
    ///     .process_pair(&pair)?
    ///     .expect("this fixture has an acceptable overlap");
    ///
    /// assert_eq!(merged.id(), "read-1");
    /// assert_eq!(merged.sequence_bytes().len(), merged.quality_bytes().len());
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// A pair with no acceptable overlap returns `Ok(None)` rather than an error.
    ///
    /// # Errors
    ///
    /// Returns an error if pairing, overlap discovery, validation, merging, or
    /// correction fail for this input pair. A successfully searched pair with no
    /// overlap returns `Ok(None)`.
    pub fn process_pair<R>(&self, pair: &PairInput<R>) -> Result<Option<OwnedSequenceRead>>
    where
        R: SeqRecordView,
    {
        self.on_pair(pair)?
            .find_overlap()?
            .and_then_found(|overlap| overlap.validate()?.merge()?.correct()?.into_owned_read())
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
        }
    }

    /// Process an iterator of paired records with a custom per-pair pipeline closure.
    ///
    /// This advanced entrypoint enables callers to choose explicit branch ordering
    /// for each pair while preserving the same per-item `Result` behavior as
    /// [`Assembler::process_iter`]. Use this when you intentionally need a non-default
    /// operation order, unvalidated intermediate slices, or mixed terminal output types.
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
    /// ```rust
    /// use libpairassembly::prelude::*;
    ///
    /// # fn main() -> libpairassembly::Result<()> {
    /// let pair = PairInput::new(
    ///     SequenceRead::try_new(
    ///         "read-1",
    ///         "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
    ///         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ///     )?,
    ///     SequenceRead::try_new(
    ///         "read-1",
    ///         "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
    ///         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ///     )?,
    /// );
    /// let assembler = Assembler::builder().build()?;
    ///
    /// let merged = assembler
    ///     .on_pair(&pair)?
    ///     .find_overlap()?
    ///     .and_then_found(|overlap| overlap.validate()?.merge()?.correct()?.into_owned_read())?;
    ///
    /// assert!(merged.is_some());
    /// # Ok(())
    /// # }
    /// ```
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
            overlap: (),
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
    pub fn with_overlap_params(mut self, overlap: OverlapParams) -> Self {
        self.config.overlap = overlap;
        self
    }

    /// Set overlap validation parameters.
    #[must_use]
    pub fn with_validator(mut self, validator: OverlapValidator) -> Self {
        self.config.validator = validator;
        self
    }

    /// Set merge-stage parameters.
    #[must_use]
    pub fn with_merge_params(mut self, merge: MergeParams) -> Self {
        self.config.merge = merge;
        self
    }

    /// Set quality-correction parameters.
    #[must_use]
    pub fn with_correction_params(mut self, correction: CorrectionParams) -> Self {
        self.config.correction = correction;
        self
    }

    /// Finalize and return the configured assembler.
    ///
    /// # Errors
    ///
    /// This currently does not fail in practice, but returns `Result` to keep
    /// construction compatible with configuration validation during build.
    pub fn build(self) -> Result<Assembler> {
        Ok(Assembler {
            config: self.config,
        })
    }
}
