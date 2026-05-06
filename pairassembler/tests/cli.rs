use std::{error::Error, fs, path::Path, process::Output};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::TempDir;

mod support;

#[test]
fn no_args_prints_long_help() -> Result<(), Box<dyn Error>> {
    let output = pairasm_command()?.output()?;

    assert_success(&output);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("PairAssembler"));
    assert!(stdout.contains("Usage: pairasm -1 <R1.fastq[.gz]> -2 <R2.fastq[.gz]> [OPTIONS]"));
    assert!(stdout.contains("--no-validate"));

    Ok(())
}

#[test]
fn mixed_run_writes_outputs_and_summary() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let pairs = support::mixed_pairs();
    let (r1, r2) = support::write_fastq_pair_files(temp.path(), "mixed", &pairs)?;
    let merged = temp.path().join("merged.fastq");
    let unmerged = temp.path().join("unmerged.fastq");
    let summary = temp.path().join("summary.json");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&merged)
        .arg("--unmerged-out")
        .arg(&unmerged)
        .arg("--summary")
        .arg(&summary)
        .arg("--progress-every")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert_success(&output);
    assert_eq!(support::count_fastq_records(&merged)?, 1);
    assert_eq!(support::count_fastq_records(&unmerged)?, 4);

    let summary_json: Value = serde_json::from_slice(&fs::read(summary)?)?;
    assert_eq!(json_u64(&summary_json, "/stats/pairs_seen")?, 3);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_merged")?, 1);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_unmerged")?, 2);
    assert_eq!(
        json_u64(&summary_json, "/stats/pairs_unmerged_no_overlap")?,
        1
    );
    assert_eq!(
        json_u64(&summary_json, "/stats/pairs_unmerged_validation_rejected")?,
        1
    );

    Ok(())
}

#[test]
fn no_validate_merges_detected_overlap_without_informativeness_check() -> Result<(), Box<dyn Error>>
{
    let temp = TempDir::new()?;
    let pairs = support::mixed_pairs();
    let (r1, r2) = support::write_fastq_pair_files(temp.path(), "no_validate", &pairs)?;
    let merged = temp.path().join("merged.fastq");
    let unmerged = temp.path().join("unmerged.fastq");
    let summary = temp.path().join("summary.json");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&merged)
        .arg("--unmerged-out")
        .arg(&unmerged)
        .arg("--summary")
        .arg(&summary)
        .arg("--no-validate")
        .arg("--progress-every")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert_success(&output);
    assert_eq!(support::count_fastq_records(&merged)?, 2);
    assert_eq!(support::count_fastq_records(&unmerged)?, 2);

    let summary_json: Value = serde_json::from_slice(&fs::read(summary)?)?;
    assert_eq!(json_u64(&summary_json, "/stats/pairs_seen")?, 3);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_merged")?, 2);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_unmerged")?, 1);
    assert_eq!(
        json_u64(&summary_json, "/stats/pairs_unmerged_validation_rejected")?,
        0
    );

    Ok(())
}

#[test]
fn gzip_inputs_and_outputs_are_supported() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let pairs = support::many_pairs(2, support::PairKind::Mergeable);
    let (r1, r2) = support::write_gzip_fastq_pair_files(temp.path(), "mergeable", &pairs)?;
    let merged = temp.path().join("merged.fastq.gz");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&merged)
        .arg("--progress-every")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert_success(&output);
    assert_eq!(support::count_fastq_records(&merged)?, 2);

    Ok(())
}

#[test]
fn merged_output_preserves_input_order() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let pairs = support::many_pairs(512, support::PairKind::Mergeable);
    let expected_names: Vec<String> = pairs.iter().map(|pair| format!("@{}", pair.id)).collect();
    let (r1, r2) = support::write_fastq_pair_files(temp.path(), "ordered", &pairs)?;
    let merged = temp.path().join("merged.fastq");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&merged)
        .arg("--progress-every")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert_success(&output);
    assert_eq!(fastq_record_names(&merged)?, expected_names);

    Ok(())
}

#[test]
fn mate_mismatch_fails_fast_without_writing_summary() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let mut pairs = support::many_pairs(1, support::PairKind::Mergeable);
    pairs[0].id = "r1-only-id".to_owned();
    let r1 = temp.path().join("mismatch_R1.fastq");
    let r2 = temp.path().join("mismatch_R2.fastq");
    support::write_fastq_pair_paths(&r1, &r2, &pairs)?;
    rewrite_first_header(&r2, "@r2-only-id/2")?;
    let summary = temp.path().join("summary.json");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("--summary")
        .arg(&summary)
        .arg("--max-mate-id-mismatches")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("paired FASTQ inputs appear to be in different orders"));
    assert!(!summary.exists());

    Ok(())
}

#[test]
fn tolerated_mate_mismatch_is_counted_and_skipped() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let pairs = support::many_pairs(2, support::PairKind::Mergeable);
    let r1 = temp.path().join("tolerated_mismatch_R1.fastq");
    let r2 = temp.path().join("tolerated_mismatch_R2.fastq");
    support::write_fastq_pair_paths(&r1, &r2, &pairs)?;
    rewrite_first_header(&r2, "@r2-only-id/2")?;
    let merged = temp.path().join("merged.fastq");
    let summary = temp.path().join("summary.json");

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&merged)
        .arg("--summary")
        .arg(&summary)
        .arg("--max-mate-id-mismatches")
        .arg("1")
        .arg("--progress-every")
        .arg("0")
        .arg("-qqq")
        .output()?;

    assert_success(&output);
    assert_eq!(support::count_fastq_records(&merged)?, 1);

    let summary_json: Value = serde_json::from_slice(&fs::read(summary)?)?;
    assert_eq!(json_u64(&summary_json, "/stats/pairs_seen")?, 2);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_processed")?, 1);
    assert_eq!(json_u64(&summary_json, "/stats/pairs_merged")?, 1);
    assert_eq!(json_u64(&summary_json, "/stats/mate_id_mismatches")?, 1);

    Ok(())
}

#[test]
fn record_count_mismatch_fails_with_input_contract_error() -> Result<(), Box<dyn Error>> {
    let temp = TempDir::new()?;
    let r1_pairs = support::many_pairs(2, support::PairKind::Mergeable);
    let r2_pairs = support::many_pairs(1, support::PairKind::Mergeable);
    let r1 = temp.path().join("count_R1.fastq");
    let r2 = temp.path().join("count_R2.fastq");
    support::write_fastq_pair_paths(&r1, &temp.path().join("discard_R2.fastq"), &r1_pairs)?;
    support::write_fastq_pair_paths(&temp.path().join("discard_R1.fastq"), &r2, &r2_pairs)?;

    let output = pairasm_command()?
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-qqq")
        .output()?;

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("paired FASTQ inputs have different record counts"));

    Ok(())
}

fn pairasm_command() -> Result<Command, Box<dyn Error>> {
    Ok(Command::cargo_bin("pairasm")?)
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "pairasm failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn json_u64(value: &Value, pointer: &str) -> Result<u64, Box<dyn Error>> {
    value
        .pointer(pointer)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing unsigned integer at JSON pointer {pointer}").into())
}

fn fastq_record_names(path: &Path) -> Result<Vec<String>, Box<dyn Error>> {
    Ok(fs::read_to_string(path)?
        .lines()
        .step_by(4)
        .map(str::to_owned)
        .collect())
}

fn rewrite_first_header(path: &Path, header: &str) -> Result<(), Box<dyn Error>> {
    let contents = fs::read_to_string(path)?;
    let mut lines = contents.lines();
    let _old_header = lines.next();
    let mut rewritten = String::new();
    rewritten.push_str(header);
    rewritten.push('\n');
    for line in lines {
        rewritten.push_str(line);
        rewritten.push('\n');
    }
    fs::write(path, rewritten)?;
    Ok(())
}
