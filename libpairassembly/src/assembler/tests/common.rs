use crate::{
    assembler::{OverlapContext, OverlapSearch, PairInput},
    test_fixtures::TupleRecord,
};

pub(super) fn rec(id: &str, seq: &str, qual: &str) -> TupleRecord {
    TupleRecord::from_strs(id, seq, qual)
}

pub(super) fn demo_pair(id: &str) -> PairInput<TupleRecord> {
    validation_demo_pair(id)
}

pub(super) fn validation_demo_pair(id: &str) -> PairInput<TupleRecord> {
    PairInput::new(
        rec(id, "ACGTTGCAGTAC", "IIIIIIIIIIII"),
        rec(id, "GTACTGCAACGT", "IIIIIIIIIIII"),
    )
}

pub(super) fn expect_found<'asm, 'pair, R>(
    search: OverlapSearch<'asm, 'pair, R>,
) -> OverlapContext<'asm, 'pair, R> {
    match search {
        OverlapSearch::Found(ctx) => ctx,
        OverlapSearch::NoOverlap(_) => panic!("fixture should have a detectable overlap"),
    }
}
