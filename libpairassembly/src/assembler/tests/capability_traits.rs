use crate::{
    assembler::{
        HasConsensusRecord, HasCorrectionEvidence, HasPairOverlap, HasReadPair, HasValidationDiag,
        OverlapContext, PairState, ValidatedContext,
    },
    merge::MergedRead,
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
    ValidatedContext<'asm, 'pair, R>: PairState + HasPairOverlap + HasReadPair + HasValidationDiag,
{
}

fn assert_validated_overlap_caps<'pair>()
where
    ValidatedOverlap<'pair>: PairState + HasPairOverlap + HasReadPair + HasValidationDiag,
{
}

fn assert_merged_read_caps()
where
    MergedRead: PairState + HasConsensusRecord + HasCorrectionEvidence,
{
}

#[test]
fn test_capability_trait_coverage_compile_assertions() {
    assert_overlap_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_context_caps::<'static, 'static, TupleRecord>();
    assert_validated_overlap_caps::<'static>();
    assert_merged_read_caps();
}
