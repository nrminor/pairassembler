use super::common::demo_pair;
use crate::assembler::{Assembler, OverlapParams};

#[test]
fn test_process_iter_with_custom_checked_merge_pipeline() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pairs = vec![demo_pair("read-custom-1"), demo_pair("read-custom-2")];

    let results = asm
        .process_iter_with(pairs, |ready| {
            let corrected = ready.overlap()?.validate()?.merge()?.correct()?;
            Ok(corrected.into_owned_read()?.sequence().to_string())
        })
        .collect::<Vec<_>>();

    assert_eq!(results.len(), 2);
    assert!(results.iter().all(Result::is_ok));
}

#[test]
fn test_process_iter_with_custom_unmerged_pipeline() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pairs = vec![demo_pair("read-custom-unmerged")];

    let result = asm
        .process_iter_with(pairs, |ready| {
            let corrected = ready.overlap()?.correct()?;
            Ok(corrected.into_owned_pair()?.id().to_string())
        })
        .next()
        .expect("iterator should yield one singleton custom-pipeline result")
        .expect("custom unvalidated pipeline should succeed for demo pair");

    assert_eq!(result, "read-custom-unmerged");
}
