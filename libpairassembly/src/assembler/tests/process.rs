use super::common::{demo_pair, rec};
use crate::{
    Error,
    assembler::{Assembler, ExecutionPolicy, OverlapParams, PairInput, TiePolicy},
    errors::OverlapError,
};

#[test]
fn test_process_pair_with_tuple_record_fixture() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read1", "ACGTACGT", "IIIIIIII"),
        rec("read1", "ACGTACGT", "IIIIIIII"),
    );

    let result = asm.process_pair(&pair);
    assert!(matches!(
        result,
        Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
    ));
}

#[test]
fn test_process_iter_yields_results() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pairs = vec![demo_pair("read1"), demo_pair("read2")];

    let results = asm.process_iter(pairs).collect::<Vec<_>>();
    assert_eq!(results.len(), 2);
}

#[test]
fn test_process_pair_equals_process_iter_singleton_success() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");

    let single = asm
        .process_pair(&demo_pair("read-single"))
        .expect("singleton process_pair should succeed for demo pair");
    let iter = asm
        .process_iter(vec![demo_pair("read-single")])
        .next()
        .expect("iterator should yield one singleton result")
        .expect("singleton process_iter result should succeed for demo pair");

    assert_eq!(single.id(), iter.id());
    assert_eq!(single.sequence_bytes(), iter.sequence_bytes());
    assert_eq!(single.quality_bytes(), iter.quality_bytes());
}

#[test]
fn test_process_pair_equals_process_iter_singleton_error() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-tie", "ACGTACGT", "IIIIIIII"),
        rec("read-tie", "ACGTACGT", "IIIIIIII"),
    );

    let single = asm.process_pair(&pair).unwrap_err();
    assert!(matches!(
        single,
        Error::OverlapError(OverlapError::OverlapTie { .. })
    ));

    let iter = asm
        .process_iter(vec![PairInput::new(
            rec("read-tie", "ACGTACGT", "IIIIIIII"),
            rec("read-tie", "ACGTACGT", "IIIIIIII"),
        )])
        .next()
        .expect("iterator should yield one singleton error result")
        .unwrap_err();
    assert!(matches!(
        iter,
        Error::OverlapError(OverlapError::OverlapTie { .. })
    ));
}

#[test]
fn test_process_iter_batch_policy_matches_record_policy() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm_record = Assembler::builder()
        .overlap(overlap)
        .execution(ExecutionPolicy::record())
        .build()
        .expect("record execution policy should build successfully");
    let asm_batch = Assembler::builder()
        .overlap(overlap)
        .execution(ExecutionPolicy::batch())
        .build()
        .expect("batch execution policy should build successfully");

    let record = asm_record
        .process_iter(vec![demo_pair("read-policy")])
        .next()
        .expect("record-policy iterator should yield a singleton result")
        .expect("record-policy singleton result should succeed");
    let batch = asm_batch
        .process_iter(vec![demo_pair("read-policy")])
        .next()
        .expect("batch-policy iterator should yield a singleton result")
        .expect("batch-policy singleton result should succeed");

    assert_eq!(record.id(), batch.id());
    assert_eq!(record.sequence_bytes(), batch.sequence_bytes());
    assert_eq!(record.quality_bytes(), batch.quality_bytes());
}

#[test]
fn test_process_pair_reports_no_overlap_outcome_at_merge_stage() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap-process", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap-process", "CCCCCCCC", "IIIIIIII"),
    );

    assert!(matches!(
        asm.process_pair(&pair),
        Err(Error::OverlapError(OverlapError::NoOverlapFound))
    ));
}

#[test]
fn test_process_iter_singleton_no_overlap_matches_process_pair_error() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
    );

    let single = asm.process_pair(&pair).unwrap_err();
    let iter = asm
        .process_iter(vec![PairInput::new(
            rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
            rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
        )])
        .next()
        .expect("iterator should yield one singleton error result")
        .unwrap_err();

    assert!(matches!(
        single,
        Error::OverlapError(OverlapError::NoOverlapFound)
    ));
    assert!(matches!(
        iter,
        Error::OverlapError(OverlapError::NoOverlapFound)
    ));
}
