use crate::{
    assembler::{
        AssemblyContext, Corrected, CorrectedContext, CorrectedMergedContext, HasConsensusRecord,
        HasPairOverlap, HasValidationMetrics, Merged, NoOverlapContext, NoOverlapFound,
        OverlapContext, OverlapFound, OverlapUnsearched, PairReady, PairState, Uncorrected,
        Unmerged, Unvalidated, ValidatedContext, ValidatedCorrectedContext,
        ValidatedCorrectedMergedContext, ValidatedMergedContext,
    },
    test_fixtures::TupleRecord,
    validate::ValidatedOverlap,
};

fn assert_overlap_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    OverlapContext<'asm, 'pair, R>: PairState
        + AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Uncorrected,
        > + HasPairOverlap,
{
}

fn assert_pair_ready_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    PairReady<'asm, 'pair, R>: AssemblyContext<
            OverlapState = OverlapUnsearched,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Uncorrected,
        >,
{
}

fn assert_validated_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    ValidatedContext<'asm, 'pair, R>: PairState + HasPairOverlap + HasValidationMetrics,
{
}

fn assert_corrected_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    CorrectedContext<'asm, 'pair, R>: PairState
        + AssemblyContext<
            OverlapState = OverlapFound,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Corrected,
        >,
{
}

fn assert_validated_corrected_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    ValidatedCorrectedContext<'asm, 'pair, R>: PairState + HasValidationMetrics,
{
}

fn assert_validated_overlap_caps<'pair>()
where
    ValidatedOverlap<'pair>: PairState + HasPairOverlap + HasValidationMetrics,
{
}

fn assert_validated_merged_context_caps<'asm, 'pair>()
where
    ValidatedMergedContext<'asm, 'pair>: PairState
        + AssemblyContext<
            OverlapState = OverlapFound,
            MergeState = Merged,
            CorrectionState = Uncorrected,
        > + HasConsensusRecord
        + HasValidationMetrics,
{
}

fn assert_corrected_merged_context_caps<'asm>()
where
    CorrectedMergedContext<'asm>: PairState
        + AssemblyContext<
            OverlapState = OverlapFound,
            MergeState = Merged,
            CorrectionState = Corrected,
        > + HasConsensusRecord,
{
}

fn assert_no_overlap_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    NoOverlapContext<'asm, 'pair, R>: PairState
        + AssemblyContext<
            OverlapState = NoOverlapFound,
            ValidationState = Unvalidated,
            MergeState = Unmerged,
            CorrectionState = Uncorrected,
        >,
{
}

fn assert_validated_corrected_merged_context_caps<'asm>()
where
    ValidatedCorrectedMergedContext<'asm>: PairState + HasConsensusRecord + HasValidationMetrics,
{
}

#[test]
fn test_capability_trait_coverage_compile_assertions() {
    assert_pair_ready_context_caps::<'static, 'static, TupleRecord>();
    assert_no_overlap_context_caps::<'static, 'static, TupleRecord>();
    assert_overlap_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_context_caps::<'static, 'static, TupleRecord>();
    assert_corrected_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_corrected_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_overlap_caps::<'static>();
    assert_validated_merged_context_caps::<'static, 'static>();
    assert_corrected_merged_context_caps::<'static>();
    assert_validated_corrected_merged_context_caps::<'static>();
}
