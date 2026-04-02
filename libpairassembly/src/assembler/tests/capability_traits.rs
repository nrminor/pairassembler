use crate::{
    assembler::{
        HasConsensusRecord, HasCorrectionWindow, HasPairOverlap, HasReadPair, HasValidationMetrics,
        OverlapContext, PairState, ValidatedContext, ValidatedMergedContext,
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

#[test]
fn test_capability_trait_coverage_compile_assertions() {
    assert_overlap_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_overlap_caps::<'static>();
    assert_validated_merged_context_caps::<'static>();
}
