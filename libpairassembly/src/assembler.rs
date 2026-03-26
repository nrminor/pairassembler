use crate::{
    BaseCallValidator, OverlapParams, ReadPair, Result, SequenceRead, TiePolicy,
    correct::{CorrectedMergedRead, CorrectionParams},
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
#[derive(Debug, Clone, Copy)]
pub enum ExecutionPolicy {
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

impl Default for ExecutionPolicy {
    fn default() -> Self {
        Self::Record
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
    /// # Errors
    ///
    /// Returns an error if pairing, overlap discovery, validation, merging, or
    /// correction fail for this input pair.
    pub fn process_pair<R>(&self, pair: PairInput<R>) -> Result<CorrectedMergedRead>
    where
        R: PairRecord,
    {
        self.on_pair(pair)?.overlap()?.validate()?.merge()
    }

    /// Process an iterator of paired records with this assembler configuration.
    ///
    /// Each output item corresponds to one input pair and is returned as a
    /// `Result` so callers can decide whether to fail-fast or handle per-pair
    /// errors inline.
    pub fn process_iter<'asm, I, R>(
        &'asm self,
        pairs: I,
    ) -> impl Iterator<Item = Result<CorrectedMergedRead>> + 'asm
    where
        I: IntoIterator<Item = PairInput<R>> + 'asm,
        R: PairRecord + 'asm,
    {
        pairs.into_iter().map(|pair| self.process_pair(pair))
    }

    /// Begin explicit per-pair processing flow.
    ///
    /// This entrypoint exists to preserve fluent per-pair APIs while collection
    /// APIs are layered on top.
    ///
    /// # Errors
    ///
    /// Returns an error only if pair initialization fails.
    pub fn on_pair<R>(&self, pair: PairInput<R>) -> Result<PairReady<'_, R>>
    where
        R: PairRecord,
    {
        Ok(PairReady {
            assembler: self,
            pair,
        })
    }
}

/// Builder for [`Assembler`].
#[derive(Debug, Clone)]
pub struct AssemblerBuilder {
    config: AssemblerConfig,
}

impl Default for AssemblerBuilder {
    fn default() -> Self {
        Self {
            config: AssemblerConfig::default(),
        }
    }
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

/// Pair-scoped processing handle returned by [`Assembler::on_pair`].
#[derive(Debug)]
pub struct PairReady<'asm, R> {
    assembler: &'asm Assembler,
    pair: PairInput<R>,
}

impl<'asm, R> PairReady<'asm, R>
where
    R: PairRecord,
{
    /// Detect overlap and enter overlap context.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap discovery fails or no overlap is found.
    pub fn overlap(self) -> Result<OverlapContext<'asm, R>> {
        let PairReady { assembler, pair } = self;
        let mates = pair.try_into_read_pair()?;
        let _overlap = mates
            .overlap(&assembler.config.overlap)?
            .ok_or_else(|| anyhow::anyhow!("No overlap found for paired reads"))?;

        Ok(OverlapContext { assembler, pair })
    }

    /// Process this pair to completion using the parent assembler.
    ///
    /// # Errors
    ///
    /// Returns any pipeline error encountered while processing this pair.
    pub fn process(self) -> Result<CorrectedMergedRead> {
        self.assembler.process_pair(self.pair)
    }
}

/// Context after overlap discovery, used as the per-pair DAG branching node.
#[derive(Debug)]
pub struct OverlapContext<'asm, R> {
    assembler: &'asm Assembler,
    pair: PairInput<R>,
}

impl<'asm, R> OverlapContext<'asm, R>
where
    R: PairRecord,
{
    /// Validate overlap with configured validator and enter validated context.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap validation fails.
    pub fn validate(self) -> Result<ValidatedContext<'asm, R>> {
        let OverlapContext { assembler, pair } = self;

        let mates = pair.try_into_read_pair()?;
        let overlap = mates
            .overlap(&assembler.config.overlap)?
            .ok_or_else(|| anyhow::anyhow!("No overlap found for paired reads"))?;
        overlap.validate(&mates, &assembler.config.validator)?;

        Ok(ValidatedContext { assembler, pair })
    }

    /// Merge and correct this overlap without validation checks.
    ///
    /// # Errors
    ///
    /// Returns an error if merge or correction fails.
    pub fn merge_unchecked(self) -> Result<CorrectedMergedRead> {
        let OverlapContext { assembler, pair } = self;

        let mates = pair.try_into_read_pair()?;
        let overlap = mates
            .overlap(&assembler.config.overlap)?
            .ok_or_else(|| anyhow::anyhow!("No overlap found for paired reads"))?;
        let validated = crate::ValidatedOverlap {
            mates: &mates,
            overlap,
        };

        let _merge_params = assembler.config.merge;
        let _correction_params = assembler.config.correction;

        validated.merge()?.correct()
    }
}

/// Context after explicit overlap validation.
#[derive(Debug)]
pub struct ValidatedContext<'asm, R> {
    assembler: &'asm Assembler,
    pair: PairInput<R>,
}

impl<R> ValidatedContext<'_, R>
where
    R: PairRecord,
{
    /// Merge and correct this pair using the checked path.
    ///
    /// # Errors
    ///
    /// Returns an error if overlap, merge, or correction fail.
    pub fn merge(self) -> Result<CorrectedMergedRead> {
        let ValidatedContext { assembler, pair } = self;

        let _merge_params = assembler.config.merge;
        let _correction_params = assembler.config.correction;

        let mates = pair.try_into_read_pair()?;
        let overlap = mates
            .overlap(&assembler.config.overlap)?
            .ok_or_else(|| anyhow::anyhow!("No overlap found for paired reads"))?;

        overlap
            .validate(&mates, &assembler.config.validator)?
            .merge()?
            .correct()
    }
}

/// Boundary trait for pair records accepted by the assembler API.
///
/// Implement this for external record types to use `Assembler` directly.
pub trait PairRecord {
    fn id(&self) -> &str;
    fn seq(&self) -> &str;
    fn qual(&self) -> &str;
}

impl PairRecord for SequenceRead<'_> {
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

impl PairRecord for (&str, &str, &str) {
    fn id(&self) -> &str {
        self.0
    }

    fn seq(&self) -> &str {
        self.1
    }

    fn qual(&self) -> &str {
        self.2
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
        R: PairRecord,
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

    #[test]
    fn test_builder_with_defaults() {
        let asm = Assembler::builder().build().unwrap();
        assert!(matches!(asm.config().execution, ExecutionPolicy::Record));
    }

    #[test]
    fn test_process_pair_with_tuple_record() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder().overlap(overlap).build().unwrap();
        let pair = PairInput::new(
            ("read1", "TTTTACGTACGT", "IIIIIIIIIIII"),
            ("read1", "ACGTACGT", "IIIIIIII"),
        );

        let result = asm.process_pair(pair);
        assert!(matches!(
            result,
            Err(crate::Error::OverlapError(
                crate::errors::OverlapError::OverlapTie(_)
            ))
        ));
    }

    #[test]
    fn test_process_iter_yields_results() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder().overlap(overlap).build().unwrap();
        let pairs = vec![
            PairInput::new(
                ("read1", "TTTTACGTACGT", "IIIIIIIIIIII"),
                ("read1", "ACGTACGT", "IIIIIIII"),
            ),
            PairInput::new(
                ("read2", "TTTTACGTACGT", "IIIIIIIIIIII"),
                ("read2", "ACGTACGT", "IIIIIIII"),
            ),
        ];

        let results = asm.process_iter(pairs).collect::<Vec<_>>();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_on_pair_process_delegates() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3)
            .with_tie_policy(TiePolicy::Reject);
        let asm = Assembler::builder().overlap(overlap).build().unwrap();
        let pair = PairInput::new(
            ("read1", "TTTTACGTACGT", "IIIIIIIIIIII"),
            ("read1", "ACGTACGT", "IIIIIIII"),
        );

        let delegated = asm.on_pair(pair).unwrap().process();
        assert!(matches!(
            delegated,
            Err(crate::Error::OverlapError(
                crate::errors::OverlapError::OverlapTie(_)
            ))
        ));
    }

    #[test]
    fn test_context_checked_and_unchecked_paths_exist() {
        let overlap = OverlapParams::default()
            .with_min_overlap(3)
            .with_min_comparisons(3);
        let asm = Assembler::builder().overlap(overlap).build().unwrap();
        let pair1 = PairInput::new(
            ("read1", "TTTTACGTACGT", "IIIIIIIIIIII"),
            ("read1", "ACGTACGT", "IIIIIIII"),
        );
        let pair2 = PairInput::new(
            ("read2", "TTTTACGTACGT", "IIIIIIIIIIII"),
            ("read2", "ACGTACGT", "IIIIIIII"),
        );

        let checked = asm.on_pair(pair1).unwrap().overlap().unwrap().validate();
        assert!(checked.is_ok());

        let unchecked = asm
            .on_pair(pair2)
            .unwrap()
            .overlap()
            .unwrap()
            .merge_unchecked();
        assert!(unchecked.is_ok());
    }
}
