use super::common::rec;
use crate::{
    Error,
    assembler::{Assembler, OverlapParams, OverlapSearch, PairInput, TiePolicy},
    errors::OverlapError,
};

#[test]
fn test_no_overlap_outcome_is_successful_search_branch() {
    let overlap = OverlapParams::default()
        .with_min_overlap(4)
        .with_min_comparisons(4);
    let asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-no-overlap", "AAAAAAAA", "IIIIIIII"),
        rec("read-no-overlap", "CCCCCCCC", "IIIIIIII"),
    );

    let search = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .find_overlap()
        .expect("overlap stage should run without scanner/conversion errors");

    let mut inspected_id = None;
    let search = search.inspect_no_overlap(|ctx| {
        inspected_id = Some(ctx.read_pair().fwd_id().to_string());
    });

    assert!(search.is_no_overlap());
    assert!(!search.is_found());
    assert_eq!(inspected_id.as_deref(), Some("read-no-overlap"));
    assert!(matches!(search, OverlapSearch::NoOverlap(_)));
}

#[test]
fn test_overlap_tie_still_errors_at_overlap_stage() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = PairInput::new(
        rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
        rec("read-tie-direct", "ACGTACGT", "IIIIIIII"),
    );

    assert!(matches!(
        asm.on_pair(&pair)
            .expect("on_pair should convert tuple records into read-pair context")
            .find_overlap(),
        Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
    ));
}
