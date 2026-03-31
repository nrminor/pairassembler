use super::common::{demo_pair, rec};
use crate::{
    Error,
    assembler::{
        Assembler, BaseCallValidator, OverlapParams, PairInput, TiePolicy, ValidatedContext,
    },
    errors::OverlapError,
};

#[test]
fn test_on_pair_process_delegates() {
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

    let delegated = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .process();
    assert!(matches!(
        delegated,
        Err(Error::OverlapError(OverlapError::OverlapTie { .. }))
    ));
}

#[test]
fn test_context_checked_and_unchecked_paths_exist() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair1 = demo_pair("read1");
    let pair2 = demo_pair("read2");

    let checked = asm
        .on_pair(&pair1)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .validate();
    assert!(checked.is_ok());

    let unchecked = asm
        .on_pair(&pair2)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .merge_unchecked();
    assert!(unchecked.is_ok());
}

#[test]
fn test_overlap_context_clone_branches_without_recomputing_overlap_selection() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = demo_pair("read-clone");

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    let checked = ctx
        .clone()
        .validate()
        .expect("validation should succeed for overlap-clone fixture")
        .merge()
        .expect("checked merge should succeed for overlap-clone fixture")
        .correct()
        .expect("checked correction should succeed for overlap-clone fixture");
    let unchecked = ctx
        .merge_unchecked()
        .expect("unchecked merge should succeed for overlap-clone fixture")
        .correct()
        .expect("unchecked correction should succeed for overlap-clone fixture");

    assert_eq!(checked.id(), unchecked.id());
    assert_eq!(checked.sequence_bytes(), unchecked.sequence_bytes());
    assert_eq!(checked.quality_bytes(), unchecked.quality_bytes());
}

#[test]
fn test_correct_pair_checked_and_unchecked_paths_match() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = demo_pair("read-correct");

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    let checked = ctx
        .clone()
        .validate()
        .expect("validation should succeed for checked-vs-unchecked fixture")
        .correct_pair()
        .expect("checked correction should succeed for checked-vs-unchecked fixture");
    let unchecked = ctx
        .correct_pair_unchecked()
        .expect("unchecked correction should succeed for checked-vs-unchecked fixture");

    assert_eq!(checked.id(), unchecked.id());
    assert_eq!(checked.fwd_sequence_bytes(), unchecked.fwd_sequence_bytes());
    assert_eq!(checked.fwd_quality_bytes(), unchecked.fwd_quality_bytes());
    assert_eq!(checked.rev_sequence_bytes(), unchecked.rev_sequence_bytes());
    assert_eq!(checked.rev_quality_bytes(), unchecked.rev_quality_bytes());
}

#[test]
fn test_correct_pair_checked_path_fails_for_low_confidence_overlap() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let validator = BaseCallValidator::new().with_min_entropy(44);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(validator)
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");
    let pair = PairInput::new(
        rec("read-low-confidence", "ACGTACGT", "IIIIIIII"),
        rec("read-low-confidence", "TCGTACGT", "IIIIIIII"),
    );

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(ctx.clone().correct_pair_unchecked().is_ok());
    assert!(
        ctx.validate()
            .and_then(ValidatedContext::correct_pair)
            .is_err()
    );
}

#[test]
fn test_validate_predicate_matches_expected_overlap_quality() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let validator = BaseCallValidator::new().with_min_entropy(44);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(validator)
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");

    let good_pair = demo_pair("read-valid-predicate");
    let good_ctx = asm
        .on_pair(&good_pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(
        good_ctx
            .is_valid()
            .expect("predicate validation should evaluate cleanly")
    );

    let low_conf_pair = PairInput::new(
        rec("read-invalid-predicate", "ACGTACGT", "IIIIIIII"),
        rec("read-invalid-predicate", "TCGTACGT", "IIIIIIII"),
    );
    let low_conf_ctx = asm
        .on_pair(&low_conf_pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(
        !low_conf_ctx
            .is_valid()
            .expect("predicate validation should evaluate cleanly")
    );
}

#[test]
fn test_validated_context_retains_validation_metrics() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = demo_pair("read-validation-metrics");

    let overlap_ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(overlap_ctx.validation_metrics_ref().is_none());

    let validated = overlap_ctx
        .validate()
        .expect("validation should succeed for retained-metrics fixture");
    let metrics = validated
        .validation_metrics_ref()
        .expect("validated contexts should retain validation metrics");

    assert!(metrics.overlap_len() >= metrics.min_overlap_len());
    assert!(metrics.mismatch_count() <= metrics.overlap_len());
}
