use crate::{
    assembler::{
        CorrectedContext, CorrectedMergedContext, HasConsensusRecord, HasCorrectionWindow,
        HasPairOverlap, HasReadPair, HasValidationMetrics, OverlapContext, PairState,
        ValidatedContext, ValidatedCorrectedContext, ValidatedCorrectedMergedContext,
        ValidatedMergedContext,
    },
    test_fixtures::TupleRecord,
    validate::ValidatedOverlap,
};

fn assert_overlap_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    OverlapContext<'asm, 'pair, R>: PairState + HasPairOverlap + HasReadPair,
{
}

fn assert_validated_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    ValidatedContext<'asm, 'pair, R>:
        PairState + HasPairOverlap + HasReadPair + HasValidationMetrics,
{
}

fn assert_corrected_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    CorrectedContext<'asm, 'pair, R>: PairState + HasPairOverlap + HasReadPair,
{
}

fn assert_validated_corrected_context_caps<'asm, 'pair, R>()
where
    R: 'pair,
    ValidatedCorrectedContext<'asm, 'pair, R>:
        PairState + HasPairOverlap + HasReadPair + HasValidationMetrics,
{
}

fn assert_validated_overlap_caps<'pair>()
where
    ValidatedOverlap<'pair>: PairState + HasPairOverlap + HasReadPair + HasValidationMetrics,
{
}

fn assert_validated_merged_context_caps<'asm>()
where
    ValidatedMergedContext<'asm>:
        PairState + HasConsensusRecord + HasCorrectionWindow + HasValidationMetrics,
{
}

fn assert_corrected_merged_context_caps<'asm>()
where
    CorrectedMergedContext<'asm>: PairState + HasConsensusRecord,
{
}

fn assert_validated_corrected_merged_context_caps<'asm>()
where
    ValidatedCorrectedMergedContext<'asm>: PairState + HasConsensusRecord + HasValidationMetrics,
{
}

#[test]
fn test_capability_trait_coverage_compile_assertions() {
    assert_overlap_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_context_caps::<'static, 'static, TupleRecord>();
    assert_corrected_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_corrected_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_overlap_caps::<'static>();
    assert_validated_merged_context_caps::<'static>();
    assert_corrected_merged_context_caps::<'static>();
    assert_validated_corrected_merged_context_caps::<'static>();
}
