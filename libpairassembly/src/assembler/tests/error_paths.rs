use super::common::rec;
use crate::{
    Error,
    assembler::{Assembler, OverlapParams, PairInput, TiePolicy},
    errors::OverlapError,
};

#[test]
fn test_no_overlap_outcome_flows_through_context_and_fails_on_consumers() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap", "CCCCCCCC", "IIIIIIII"),
    );

    let overlapped = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    let validated = overlapped
        .clone()
        .validate()
        .expect("validation should succeed even when no-overlap is carried forward");

    assert!(matches!(
        overlapped.clone().merge_unchecked(),
        Err(Error::OverlapError(OverlapError::NoOverlapFound))
    ));
    assert!(matches!(
        validated.clone().merge(),
        Err(Error::OverlapError(OverlapError::NoOverlapFound))
    ));
    assert!(matches!(
        overlapped.correct_pair_unchecked(),
        Err(Error::OverlapError(OverlapError::NoOverlapFound))
    ));
    assert!(matches!(
        validated.correct_pair(),
        Err(Error::OverlapError(OverlapError::NoOverlapFound))
    ));
}

#[test]
fn test_overlap_tie_still_errors_at_overlap_stage() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
        rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
    );

    assert!(matches!(
        asm.on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .overlap(),
        Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
    ));
}
