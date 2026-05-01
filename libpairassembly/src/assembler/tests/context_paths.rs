use super::common::{demo_pair, rec, validation_demo_pair};
use crate::{
    Error,
    assembler::{
        Assembler, BaseCallValidator, HasValidationMetrics, OverlapParams, PairInput, TiePolicy,
        ValidatedContext,
    },
    errors::OverlapError,
    prelude::utils::reverse_complement,
    validate::ValidationPreset,
};

#[test]
fn test_on_pair_process_delegates() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3)
        .with_tie_policy(TiePolicy::Reject);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(BaseCallValidator::from_preset(ValidationPreset::Loose))
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

fn relaxed_loose_validator() -> BaseCallValidator {
    BaseCallValidator::from_preset(ValidationPreset::Loose)
        .with_k(1)
        .with_min_complexity_score(4)
}

#[test]
fn test_borrowed_read_egress_uses_fastq_ascii_qualities() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-egress");

    let ready = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context");
    let ready_pair = ready.as_read_pair();
    assert!(
        ready_pair
            .fwd_quality_scores()
            .bytes()
            .all(|quality| quality == b'I')
    );

    let uncorrected_merged = ready
        .clone()
        .overlap()
        .expect("overlap stage should run before uncorrected merge egress")
        .merge()
        .expect("merge should succeed for borrowed uncorrected egress fixture");
    let uncorrected_merged_read = uncorrected_merged.as_merged_read();
    assert_eq!(
        uncorrected_merged_read.sequence().len(),
        uncorrected_merged_read.quality_scores().len()
    );
    assert!(
        uncorrected_merged_read
            .quality_scores()
            .bytes()
            .all(|quality| quality >= b'!')
    );
}

#[test]
fn test_chainable_owned_egress_hides_context_carriers() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");

    let merged = asm
        .on_pair(&validation_demo_pair("read-owned-merged"))
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .validate()
        .expect("validation should succeed before owned merged egress")
        .merge()
        .expect("merge should succeed before owned merged egress")
        .correct()
        .expect("correction should succeed before owned merged egress")
        .into_owned_read()
        .expect("corrected merged context should convert into an owned read");
    assert_eq!(merged.id(), "read-owned-merged");
    assert_eq!(merged.sequence().len(), merged.quality_scores().len());
    assert!(
        merged
            .quality_scores()
            .bytes()
            .all(|quality| quality >= b'!')
    );

    let pair = asm
        .on_pair(&validation_demo_pair("read-owned-pair"))
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .correct()
        .expect("correction should succeed before owned pair egress")
        .into_owned_pair()
        .expect("corrected pair context should convert into an owned pair");
    let borrowed = pair.as_read_pair();
    assert_eq!(borrowed.fwd_id(), "read-owned-pair");
    assert_eq!(borrowed.rev_id(), "read-owned-pair");
    assert_eq!(
        borrowed.fwd_sequence().len(),
        borrowed.fwd_quality_scores().len()
    );
    assert_eq!(
        borrowed.rev_sequence().len(),
        borrowed.rev_quality_scores().len()
    );
}

#[test]
fn test_context_validated_and_unvalidated_paths_exist() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair1 = validation_demo_pair("read1");
    let pair2 = validation_demo_pair("read2");

    let checked = asm
        .on_pair(&pair1)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .validate();
    assert!(checked.is_ok());

    let unvalidated = asm
        .on_pair(&pair2)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .merge();
    assert!(unvalidated.is_ok());
}

#[test]
fn test_overlap_context_clone_branches_without_recomputing_overlap_selection() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-clone");

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
    let unvalidated = ctx
        .merge()
        .expect("unvalidated merge should succeed for overlap-clone fixture")
        .correct()
        .expect("merged correction should succeed for overlap-clone fixture");
    let checked = checked
        .into_owned_read()
        .expect("checked corrected merge should convert to owned read");
    let unvalidated = unvalidated
        .into_owned_read()
        .expect("unvalidated corrected merge should convert to owned read");

    assert_eq!(checked.id(), unvalidated.id());
    assert_eq!(checked.sequence(), unvalidated.sequence());
    assert_eq!(checked.quality_scores(), unvalidated.quality_scores());
}

#[test]
fn test_correct_pair_validated_and_unvalidated_paths_match() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-correct");

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    let checked = ctx
        .clone()
        .validate()
        .expect("validation should succeed for validated-vs-unvalidated fixture")
        .correct()
        .expect("validated correction should succeed for validated-vs-unvalidated fixture");
    let unvalidated = ctx
        .correct()
        .expect("unvalidated correction should succeed for validated-vs-unvalidated fixture");
    let checked = checked
        .into_owned_pair()
        .expect("validated corrected pair should convert to owned pair");
    let unvalidated = unvalidated
        .into_owned_pair()
        .expect("unvalidated corrected pair should convert to owned pair");
    let checked = checked.as_read_pair();
    let unvalidated = unvalidated.as_read_pair();

    assert_eq!(checked.fwd_id(), unvalidated.fwd_id());
    assert_eq!(checked.fwd_sequence(), unvalidated.fwd_sequence());
    assert_eq!(checked.fwd_quality_bytes(), unvalidated.fwd_quality_bytes());
    assert_eq!(checked.rev_sequence(), unvalidated.rev_sequence());
    assert_eq!(checked.rev_quality_bytes(), unvalidated.rev_quality_bytes());
}

#[test]
fn test_unvalidated_pair_correction_keeps_overlap_reverse_complement_consistent() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-correct-overlap-consistency");

    let corrected = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .correct()
        .expect("unvalidated pair correction should succeed for correction-consistency fixture");
    let corrected = corrected
        .into_owned_pair()
        .expect("corrected pair should convert to owned pair");
    let corrected = corrected.as_read_pair();

    let rev_rc = reverse_complement(corrected.rev_sequence());
    assert_eq!(corrected.fwd_sequence(), rev_rc);
}

#[test]
fn test_correct_pair_checked_path_fails_for_low_confidence_overlap() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let validator = BaseCallValidator::from_preset(ValidationPreset::Strict);
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
    assert!(ctx.clone().correct().is_ok());
    assert!(ctx.validate().and_then(ValidatedContext::correct).is_err());
}

#[test]
fn test_corrected_pair_context_validates_corrected_evidence() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let validator = BaseCallValidator::from_preset(ValidationPreset::Normal)
        .with_k(1)
        .with_min_complexity_score(4);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(validator)
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");
    let pair = PairInput::new(
        rec("read-correct-then-validate", "ACGTACGT", "IIIIIIII"),
        rec("read-correct-then-validate", "TCGTACGT", "IIIIIIII"),
    );

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");

    assert!(ctx.clone().validate().is_err());

    let validated_corrected = ctx
        .correct()
        .expect("unvalidated correction should succeed before corrected validation")
        .validate()
        .expect("corrected evidence should validate after pair correction");
    assert_eq!(validated_corrected.validation_metrics().mismatch_count(), 0);
}

#[test]
fn test_corrected_pair_context_merges_corrected_evidence() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let validator = BaseCallValidator::from_preset(ValidationPreset::Normal)
        .with_k(1)
        .with_min_complexity_score(4);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(validator)
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");
    let pair = PairInput::new(
        rec("read-correct-merge", "ACGTACGT", "IIIIIIII"),
        rec("read-correct-merge", "TCGTACGT", "IIIIIIII"),
    );

    let ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");

    let checked = ctx
        .clone()
        .correct()
        .expect("pair correction should succeed before corrected checked merge")
        .validate()
        .expect("corrected pair evidence should validate before checked merge")
        .merge()
        .expect("checked merge should accept corrected pair evidence");
    let unvalidated = ctx
        .correct()
        .expect("pair correction should succeed before corrected unvalidated merge")
        .merge()
        .expect("unvalidated merge should accept corrected pair evidence");
    let checked = checked
        .into_owned_read()
        .expect("checked corrected merge should convert to owned read");
    let unvalidated = unvalidated
        .into_owned_read()
        .expect("unvalidated corrected merge should convert to owned read");

    assert_eq!(checked.id(), unvalidated.id());
    assert_eq!(checked.sequence(), unvalidated.sequence());
    assert_eq!(checked.quality_scores(), unvalidated.quality_scores());
}

#[test]
fn test_validate_predicate_matches_expected_overlap_quality() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let good_asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");

    let good_pair = validation_demo_pair("read-valid-predicate");
    let good_ctx = good_asm
        .on_pair(&good_pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(
        good_ctx
            .is_valid()
            .expect("predicate validation should evaluate cleanly")
    );

    let low_conf_asm = Assembler::builder()
        .overlap(overlap)
        .validate(BaseCallValidator::from_preset(ValidationPreset::Strict))
        .build()
        .expect("assembler builder should accept explicit overlap/validation settings");

    let low_conf_pair = PairInput::new(
        rec("read-invalid-predicate", "ACGTACGT", "IIIIIIII"),
        rec("read-invalid-predicate", "TCGTACGT", "IIIIIIII"),
    );
    let low_conf_ctx = low_conf_asm
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
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-validation-metrics");

    let overlap_ctx = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors");
    assert!(overlap_ctx.validation_metrics_ref().is_none());

    let validated = overlap_ctx
        .validate()
        .expect("validation should succeed for retained-metrics fixture");
    let metrics = validated.validation_metrics();

    assert!(metrics.overlap_len() >= metrics.min_informative_overlap_len());
    assert!(metrics.mismatch_count() <= metrics.overlap_len());
}

#[test]
fn test_validated_context_predicate_short_circuits_from_retained_metrics() {
    let overlap = OverlapParams::default()
        .with_min_overlap(3)
        .with_min_comparisons(3);
    let asm = Assembler::builder()
        .overlap(overlap)
        .validate(relaxed_loose_validator())
        .build()
        .expect("assembler builder should accept explicit overlap settings");
    let pair = validation_demo_pair("read-valid-short-circuit");

    let validated = asm
        .on_pair(&pair)
        .expect("on_pair should convert tuple records into read-pair context")
        .overlap()
        .expect("overlap stage should run without scanner/conversion errors")
        .validate()
        .expect("validation should succeed for short-circuit fixture");

    assert!(
        validated
            .is_valid()
            .expect("validated predicate should use retained metrics cleanly")
    );
}
