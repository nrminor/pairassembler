use crate::{assembler::PairInput, test_fixtures::TupleRecord};

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
