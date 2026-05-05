use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    process::Command,
};

use color_eyre::eyre::{Result, bail};

use crate::{
    cli::RunOptions,
    commands::build_tool_command,
    config::{read_datasets, read_subset_metadata, write_tool_versions},
    fastq::{fastq_record_count, file_size},
    model::{HyperfineReport, HyperfineResult, SubsetMetadata, Tool, ToolCommand, ToolPaths},
    process::run_command,
    shell::shell_join,
    validate::validate_tool_run,
};

pub fn run_matrix(options: &RunOptions) -> Result<()> {
    let paths = ToolPaths::from_environment()?;
    let datasets = read_datasets(&options.common.config)?;
    let run_id = utc_run_id()?;
    let run_dir = options.runs_root.join(&run_id);
    fs::create_dir_all(run_dir.join("metadata"))?;
    write_run_metadata(&run_dir, &run_id, options)?;
    write_tool_versions(&run_dir, &paths)?;

    for dataset in datasets {
        let subset =
            read_subset_metadata(&options.common.data_root, &dataset.name, options.read_pairs)?;
        for tool in &options.tools {
            run_tool(&paths, options, &run_id, &run_dir, &subset, *tool)?;
        }
    }

    eprintln!("Run artifacts: {}", run_dir.display());
    Ok(())
}

fn run_tool(
    paths: &ToolPaths,
    options: &RunOptions,
    run_id: &str,
    run_dir: &Path,
    subset: &SubsetMetadata,
    tool: Tool,
) -> Result<()> {
    let out_dir = run_dir
        .join(&subset.name)
        .join(format!("{}_pairs", subset.read_pairs))
        .join(tool.name());
    fs::create_dir_all(&out_dir)?;
    let command = build_tool_command(paths, options, subset, tool, &out_dir);
    let stdout_log = out_dir.join(format!("{}.stdout.log", tool.name()));
    let stderr_log = out_dir.join(format!("{}.stderr.log", tool.name()));
    let command_string = format!(
        "{} > {} 2> {}",
        shell_join(&command.args),
        shell_join(&[stdout_log.to_string_lossy().into_owned()]),
        shell_join(&[stderr_log.to_string_lossy().into_owned()])
    );
    fs::write(out_dir.join("command.sh"), format!("{command_string}\n"))?;

    eprintln!("[{run_id}] {} {}", subset.name, tool.name());
    let hyperfine_json = out_dir.join("hyperfine.json");
    let hyperfine_md = out_dir.join("hyperfine.md");
    run_command(
        Command::new(&paths.hyperfine)
            .arg("--runs")
            .arg(options.replicates.to_string())
            .arg("--warmup")
            .arg("1")
            .arg("--export-json")
            .arg(&hyperfine_json)
            .arg("--export-markdown")
            .arg(&hyperfine_md)
            .arg("--command-name")
            .arg(tool.name())
            .arg(&command_string),
    )?;

    let report: HyperfineReport = serde_json::from_reader(File::open(&hyperfine_json)?)?;
    let result = report
        .results
        .first()
        .ok_or_else(|| color_eyre::eyre::eyre!("hyperfine report had no results"))?;
    let merged_reads = command
        .merged_output
        .exists()
        .then(|| fastq_record_count(&command.merged_output))
        .transpose()?
        .unwrap_or(0);
    validate_tool_run(
        command.tool,
        options.mode,
        subset,
        &out_dir,
        &command.merged_output,
        merged_reads,
    )?;
    write_tool_result(
        &out_dir,
        run_id,
        options,
        subset,
        &command,
        result,
        merged_reads,
    )
}

fn write_tool_result(
    out_dir: &Path,
    run_id: &str,
    options: &RunOptions,
    subset: &SubsetMetadata,
    command: &ToolCommand,
    result: &HyperfineResult,
    merged_reads: usize,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out_dir.join("result.tsv"))?);
    writeln!(
        writer,
        "run_id\tbenchmark_mode\tdataset\taccession\tread_pairs\ttool\treplicates\tthreads\toutput_compression\tmean_s\tmedian_s\tstddev_s\tmin_s\tmax_s\tuser_s\tsystem_s\tmerged_reads\tr1_bytes\tr2_bytes\toutput_dir"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        run_id,
        options.mode,
        subset.name,
        subset.accession,
        subset.read_pairs,
        command.tool.name(),
        options.replicates,
        options.threads,
        options.output_compression,
        result.mean,
        result.median,
        optional_f64(result.stddev),
        result.min,
        result.max,
        result.user,
        result.system,
        merged_reads,
        file_size(&subset.r1)?,
        file_size(&subset.r2)?,
        out_dir.display()
    )?;
    writer.flush()?;
    Ok(())
}

fn write_run_metadata(run_dir: &Path, run_id: &str, options: &RunOptions) -> Result<()> {
    let mut writer = BufWriter::new(File::create(run_dir.join("metadata").join("run.tsv"))?);
    writeln!(writer, "key\tvalue")?;
    writeln!(writer, "run_id\t{run_id}")?;
    writeln!(writer, "benchmark_mode\t{}", options.mode)?;
    writeln!(writer, "read_pairs\t{}", options.read_pairs)?;
    writeln!(writer, "replicates\t{}", options.replicates)?;
    writeln!(writer, "threads\t{}", options.threads)?;
    writeln!(writer, "output_compression\t{}", options.output_compression)?;
    writeln!(writer, "config\t{}", options.common.config.display())?;
    writer.flush()?;
    Ok(())
}

fn optional_f64(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn utc_run_id() -> Result<String> {
    let output = Command::new("date")
        .args(["-u", "+%Y%m%dT%H%M%SZ"])
        .output()?;
    if !output.status.success() {
        bail!("failed to create UTC run id with date command");
    }
    Ok(String::from_utf8(output.stdout)?.trim().to_owned())
}
