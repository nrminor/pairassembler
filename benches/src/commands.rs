use std::path::{Path, PathBuf};

use crate::{
    cli::{BenchmarkMode, OutputCompression, RunOptions},
    model::{SubsetMetadata, Tool, ToolCommand, ToolPaths},
};

pub fn build_tool_command(
    paths: &ToolPaths,
    options: &RunOptions,
    subset: &SubsetMetadata,
    tool: Tool,
    out_dir: &Path,
) -> ToolCommand {
    let merged_output = out_dir.join(merged_output_name(tool, options.output_compression));
    let r1 = subset.r1.to_string_lossy().into_owned();
    let r2 = subset.r2.to_string_lossy().into_owned();
    let merged = merged_output.to_string_lossy().into_owned();
    let binary = paths.path_for(tool).to_string_lossy().into_owned();

    let args = match tool {
        Tool::Pairasm => pairasm_command(binary, r1, r2, merged, out_dir),
        Tool::Fastp => fastp_command(binary, r1, r2, merged, out_dir, options.threads),
        Tool::Bbmerge => bbmerge_command(binary, r1, r2, merged, out_dir, options.threads),
        Tool::Vsearch => vsearch_command(binary, r1, r2, merged, out_dir, options),
    };

    ToolCommand {
        tool,
        args,
        merged_output,
    }
}

fn pairasm_command(
    binary: String,
    r1: String,
    r2: String,
    merged: String,
    out_dir: &Path,
) -> Vec<String> {
    vec![
        binary,
        "-1".to_owned(),
        r1,
        "-2".to_owned(),
        r2,
        "-o".to_owned(),
        merged,
        "--unmerged-out".to_owned(),
        path_string(out_dir.join("pairasm.unmerged.fastq")),
        "--summary".to_owned(),
        path_string(out_dir.join("pairasm.summary.json")),
        "--progress-every".to_owned(),
        "0".to_owned(),
        "-qqq".to_owned(),
    ]
}

fn fastp_command(
    binary: String,
    r1: String,
    r2: String,
    merged: String,
    out_dir: &Path,
    threads: usize,
) -> Vec<String> {
    vec![
        binary,
        "-i".to_owned(),
        r1,
        "-I".to_owned(),
        r2,
        "--merge".to_owned(),
        "--merged_out".to_owned(),
        merged,
        "--unpaired1".to_owned(),
        path_string(out_dir.join("fastp.unpaired1.fastq")),
        "--unpaired2".to_owned(),
        path_string(out_dir.join("fastp.unpaired2.fastq")),
        "--failed_out".to_owned(),
        path_string(out_dir.join("fastp.failed.fastq")),
        "--thread".to_owned(),
        threads.to_string(),
        "--html".to_owned(),
        path_string(out_dir.join("fastp.html")),
        "--json".to_owned(),
        path_string(out_dir.join("fastp.json")),
    ]
}

fn bbmerge_command(
    binary: String,
    r1: String,
    r2: String,
    merged: String,
    out_dir: &Path,
    threads: usize,
) -> Vec<String> {
    vec![
        binary,
        "-da".to_owned(),
        format!("in1={r1}"),
        format!("in2={r2}"),
        format!("out={merged}"),
        format!(
            "outu1={}",
            out_dir.join("bbmerge.unmerged1.fastq").display()
        ),
        format!(
            "outu2={}",
            out_dir.join("bbmerge.unmerged2.fastq").display()
        ),
        format!("threads={threads}"),
        "ow=t".to_owned(),
    ]
}

fn vsearch_command(
    binary: String,
    r1: String,
    r2: String,
    merged: String,
    out_dir: &Path,
    options: &RunOptions,
) -> Vec<String> {
    let mut args = vec![
        binary,
        "--fastq_mergepairs".to_owned(),
        r1,
        "--reverse".to_owned(),
        r2,
    ];

    if options.mode == BenchmarkMode::TunedComparability {
        args.extend([
            "--fastq_allowmergestagger".to_owned(),
            "--fastq_minovlen".to_owned(),
            "30".to_owned(),
            "--fastq_maxdiffs".to_owned(),
            "5".to_owned(),
            "--fastq_maxdiffpct".to_owned(),
            "20".to_owned(),
        ]);
    }

    args.extend([
        "--fastqout".to_owned(),
        merged,
        "--fastqout_notmerged_fwd".to_owned(),
        path_string(out_dir.join("vsearch.unmerged1.fastq")),
        "--fastqout_notmerged_rev".to_owned(),
        path_string(out_dir.join("vsearch.unmerged2.fastq")),
        "--threads".to_owned(),
        options.threads.to_string(),
    ]);

    if options.mode == BenchmarkMode::TunedComparability {
        args.push("--no_progress".to_owned());
    }

    args
}

fn merged_output_name(tool: Tool, output_compression: OutputCompression) -> String {
    match output_compression {
        OutputCompression::Plain => format!("{}.merged.fastq", tool.name()),
        OutputCompression::Gzip => format!("{}.merged.fastq.gz", tool.name()),
    }
}

fn path_string(path: PathBuf) -> String {
    path.to_string_lossy().into_owned()
}
