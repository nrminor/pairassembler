use crate::{assembler::PairInput, test_fixtures::TupleRecord};

pub(super) fn rec(id: &str, seq: &str, qual: &str) -> TupleRecord {
    TupleRecord::from_strs(id, seq, qual)
}

pub(super) fn demo_pair(id: &str) -> PairInput<TupleRecord> {
    PairInput::new(rec(id, "TTTACGTA", "IIIIIIII"), rec(id, "TACGT", "IIIII"))
}
