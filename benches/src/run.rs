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
    db::{BenchmarkDb, RunRecord, ToolExecutionRecord, collect_vcs_metadata},
    fastq::{fastq_record_count, file_size},
    model::{HyperfineReport, HyperfineResult, SubsetMetadata, Tool, ToolCommand, ToolPaths},
    process::run_command,
    products::read_merged_products,
    shell::shell_join,
    validate::validate_tool_run,
};

pub fn run_matrix(options: &RunOptions) -> Result<()> {
    BenchmarkRun::new(options)?.run()
}

struct BenchmarkRun<'options> {
    options: &'options RunOptions,
    paths: ToolPaths,
    database: BenchmarkDb,
    run_key: String,
    run_id: String,
    run_dir: std::path::PathBuf,
}

impl<'options> BenchmarkRun<'options> {
    fn new(options: &'options RunOptions) -> Result<Self> {
        let paths = ToolPaths::from_environment()?;
        let run_id = utc_run_id()?;
        let run_key = uuid::Uuid::new_v4().to_string();
        let run_dir = options.runs_root.join(&run_id);
        let database = BenchmarkDb::open(&options.db)?;

        Ok(Self {
            options,
            paths,
            database,
            run_key,
            run_id,
            run_dir,
        })
    }

    fn run(self) -> Result<()> {
        let datasets = read_datasets(&self.options.common.config)?;
        fs::create_dir_all(self.run_dir.join("metadata"))?;
        write_run_metadata(&self.run_dir, &self.run_id, self.options)?;
        write_tool_versions(&self.run_dir, &self.paths)?;
        let vcs = collect_vcs_metadata();
        self.database.insert_run(&RunRecord {
            run_key: &self.run_key,
            run_label: &self.run_id,
            created_at: &self.run_id,
            run_dir: &self.run_dir,
            options: self.options,
            vcs: &vcs,
        })?;
        self.database
            .insert_tool_versions(&self.run_key, &self.paths)?;

        for dataset in datasets {
            let subset = read_subset_metadata(
                &self.options.common.data_root,
                &dataset.name,
                self.options.read_pairs,
            )?;
            self.database
                .insert_dataset_subset(&self.run_key, &subset)?;
            for tool in &self.options.tools {
                self.run_tool(&subset, *tool)?;
            }
        }

        eprintln!("Run artifacts: {}", self.run_dir.display());
        eprintln!("Benchmark database: {}", self.options.db.display());
        Ok(())
    }

    fn run_tool(&self, subset: &SubsetMetadata, tool: Tool) -> Result<()> {
        let artifacts = ToolRunArtifacts::new(&self.run_dir, subset, tool);
        fs::create_dir_all(&artifacts.out_dir)?;
        let command =
            build_tool_command(&self.paths, self.options, subset, tool, &artifacts.out_dir);
        let command_string = artifacts.command_string(&command);
        fs::write(artifacts.command_script(), format!("{command_string}\n"))?;

        eprintln!("[{}] {} {}", self.run_id, subset.name, tool.name());
        run_command(
            Command::new(&self.paths.hyperfine)
                .arg("--runs")
                .arg(self.options.replicates.to_string())
                .arg("--warmup")
                .arg("1")
                .arg("--export-json")
                .arg(&artifacts.hyperfine_json)
                .arg("--export-markdown")
                .arg(&artifacts.hyperfine_md)
                .arg("--command-name")
                .arg(tool.name())
                .arg(&command_string),
        )?;

        let report: HyperfineReport =
            serde_json::from_reader(File::open(&artifacts.hyperfine_json)?)?;
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
            self.options.mode,
            subset,
            &artifacts.out_dir,
            &command.merged_output,
            merged_reads,
        )?;
        let merged_products = read_merged_products(&command.merged_output)?;
        self.database
            .insert_merged_products(&self.run_key, subset, &command, &merged_products)?;
        let merged_product_rows = merged_products.len();
        self.database.insert_tool_execution(&ToolExecutionRecord {
            run_key: &self.run_key,
            subset,
            command: &command,
            command_string: &command_string,
            result,
            merged_reads,
            merged_product_rows,
            output_dir: &artifacts.out_dir,
        })?;
        write_tool_result(ToolResultRecord {
            out_dir: &artifacts.out_dir,
            run_id: &self.run_id,
            options: self.options,
            subset,
            command: &command,
            result,
            merged_reads,
            merged_product_rows,
        })?;
        artifacts.record_in_database(&self.database, &self.run_key, &subset.name, &command)
    }
}

struct ToolRunArtifacts {
    tool: Tool,
    out_dir: std::path::PathBuf,
    stdout_log: std::path::PathBuf,
    stderr_log: std::path::PathBuf,
    hyperfine_json: std::path::PathBuf,
    hyperfine_md: std::path::PathBuf,
}

impl ToolRunArtifacts {
    fn new(run_dir: &Path, subset: &SubsetMetadata, tool: Tool) -> Self {
        let out_dir = run_dir
            .join(&subset.name)
            .join(format!("{}_pairs", subset.read_pairs))
            .join(tool.name());
        Self {
            stdout_log: out_dir.join(format!("{}.stdout.log", tool.name())),
            stderr_log: out_dir.join(format!("{}.stderr.log", tool.name())),
            hyperfine_json: out_dir.join("hyperfine.json"),
            hyperfine_md: out_dir.join("hyperfine.md"),
            tool,
            out_dir,
        }
    }

    fn command_script(&self) -> std::path::PathBuf {
        self.out_dir.join("command.sh")
    }

    fn command_string(&self, command: &ToolCommand) -> String {
        format!(
            "{} > {} 2> {}",
            shell_join(&command.args),
            shell_join(&[self.stdout_log.to_string_lossy().into_owned()]),
            shell_join(&[self.stderr_log.to_string_lossy().into_owned()])
        )
    }

    fn record_in_database(
        &self,
        database: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        command: &ToolCommand,
    ) -> Result<()> {
        for artifact in self.artifacts(command) {
            database.insert_artifact(
                run_key,
                dataset_name,
                self.tool.name(),
                artifact.kind,
                &artifact.path,
            )?;
        }
        Ok(())
    }

    fn artifacts(&self, command: &ToolCommand) -> Vec<ToolArtifact> {
        vec![
            ToolArtifact::new("command", self.command_script()),
            ToolArtifact::new("stdout_log", self.stdout_log.clone()),
            ToolArtifact::new("stderr_log", self.stderr_log.clone()),
            ToolArtifact::new("hyperfine_json", self.hyperfine_json.clone()),
            ToolArtifact::new("hyperfine_markdown", self.hyperfine_md.clone()),
            ToolArtifact::new("merged_fastq", command.merged_output.clone()),
            ToolArtifact::new("result_tsv", self.out_dir.join("result.tsv")),
        ]
    }
}

struct ToolArtifact {
    kind: &'static str,
    path: std::path::PathBuf,
}

impl ToolArtifact {
    fn new(kind: &'static str, path: std::path::PathBuf) -> Self {
        Self { kind, path }
    }
}

struct ToolResultRecord<'a> {
    out_dir: &'a Path,
    run_id: &'a str,
    options: &'a RunOptions,
    subset: &'a SubsetMetadata,
    command: &'a ToolCommand,
    result: &'a HyperfineResult,
    merged_reads: usize,
    merged_product_rows: usize,
}

fn write_tool_result(record: ToolResultRecord<'_>) -> Result<()> {
    let mut writer = BufWriter::new(File::create(record.out_dir.join("result.tsv"))?);
    writeln!(
        writer,
        "run_id\tbenchmark_mode\tdataset\taccession\tread_pairs\ttool\treplicates\tthreads\toutput_compression\tmean_s\tmedian_s\tstddev_s\tmin_s\tmax_s\tuser_s\tsystem_s\tmerged_reads\tmerged_product_rows\tr1_bytes\tr2_bytes\toutput_dir"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        record.run_id,
        record.options.mode,
        record.subset.name,
        record.subset.accession,
        record.subset.read_pairs,
        record.command.tool.name(),
        record.options.replicates,
        record.options.threads,
        record.options.output_compression,
        record.result.mean,
        record.result.median,
        optional_f64(record.result.stddev),
        record.result.min,
        record.result.max,
        record.result.user,
        record.result.system,
        record.merged_reads,
        record.merged_product_rows,
        file_size(&record.subset.r1)?,
        file_size(&record.subset.r2)?,
        record.out_dir.display()
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
