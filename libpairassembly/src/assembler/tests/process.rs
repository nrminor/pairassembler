use super::common::{demo_pair, rec};
use crate::{
    Error,
    assembler::{Assembler, OverlapParams, PairInput, TiePolicy},
    errors::OverlapError,
};

#[test]
fn test_process_pair_with_tuple_record_fixture() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .with_overlap_params(overlap)
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
        .with_overlap_params(overlap)
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
        .with_overlap_params(overlap)
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

    let single = single.expect("demo pair should produce a merged read");
    let iter = iter.expect("demo pair should produce a merged read");

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
        .with_overlap_params(overlap)
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
fn test_process_pair_reports_no_overlap_as_empty_success() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap-process", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap-process", "CCCCCCCC", "IIIIIIII"),
    );

    assert!(
        asm.process_pair(&pair)
            .expect("no-overlap should not be an operational failure")
            .is_none()
    );
}

#[test]
fn test_process_iter_singleton_no_overlap_matches_process_pair_empty_success() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
    );

    let single = asm
        .process_pair(&pair)
        .expect("no-overlap should not be an operational failure");
    let iter = asm
        .process_iter(vec![PairInput::new(
            rec("read-no-overlap-iter", "AAAAAAAA", "IIIIIIII"),
            rec("read-no-overlap-iter", "CCCCCCCC", "IIIIIIII"),
        )])
        .next()
        .expect("iterator should yield one singleton result")
        .expect("no-overlap should not be an operational failure");

    assert!(single.is_none());
    assert!(iter.is_none());
}
