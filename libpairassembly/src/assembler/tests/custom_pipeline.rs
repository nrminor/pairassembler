use super::common::demo_pair;
use crate::assembler::{Assembler, OverlapParams};

#[test]
fn test_process_iter_with_custom_checked_merge_pipeline() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let mut asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pairs = vec![demo_pair("read-custom-1"), demo_pair("read-custom-2")];

    let results = asm
        .process_iter_with(pairs, |ready| {
            ready.find_overlap()?.and_then_found(|found| {
                let corrected = found.validate()?.merge()?.correct()?;
                Ok(corrected.into_owned_read()?.sequence().to_string())
            })
        })
        .collect::<Vec<_>>();

    assert_eq!(results.len(), 2);
    assert!(results.iter().all(|result| matches!(result, Ok(Some(_)))));
}

#[test]
fn test_process_iter_with_custom_unmerged_pipeline() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let mut asm = Assembler::builder()
        .with_overlap_params(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pairs = vec![demo_pair("read-custom-unmerged")];

    let result = asm
        .process_iter_with(pairs, |ready| {
            ready.find_overlap()?.and_then_found(|found| {
                let corrected = found.correct()?;
                Ok(corrected.into_owned_pair()?.id().to_string())
            })
        })
        .next()
        .expect("iterator should yield one singleton custom-pipeline result")
        .expect("custom unvalidated pipeline should succeed for demo pair")
        .expect("demo pair should have overlap");

    assert_eq!(result, "read-custom-unmerged");
}
