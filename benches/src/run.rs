use std::{
    fs::{self, File},
    path::Path,
    process::Command,
};

use color_eyre::eyre::{Result, WrapErr, bail};

use crate::{
    artifacts::{ArtifactKind, ArtifactRecord, ArtifactRequirement},
    cli::RunOptions,
    commands::build_tool_command,
    config::{effective_read_pairs, read_datasets, read_subset_metadata},
    db::{BenchmarkDb, RunRecord, ToolExecutionRecord},
    fastq::fastq_record_count,
    model::{HyperfineReport, SubsetMetadata, Tool, ToolCommand, ToolPaths},
    process::run_command,
    shell::shell_join,
    ui,
    validate::validate_tool_run,
    vcs::collect_vcs_metadata,
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
        let vcs = collect_vcs_metadata();
        self.database.insert_run(&RunRecord {
            run_key: &self.run_key,
            created_at: &self.run_id,
            run_dir: &self.run_dir,
            options: self.options,
            vcs: &vcs,
        })?;
        self.database
            .insert_tool_versions(&self.run_key, &self.paths)?;

        let mut execution_order = 0;
        for dataset in datasets {
            let read_pairs = effective_read_pairs(&dataset, self.options.read_pairs);
            let subset =
                read_subset_metadata(&self.options.common.data_root, &dataset.name, read_pairs)?;
            self.database
                .insert_dataset_subset(&self.run_key, &subset)?;
            for tool in &self.options.tools {
                self.run_tool(&subset, *tool, execution_order)?;
                execution_order += 1;
            }
        }

        self.database
            .mark_run_completed(&self.run_key, &utc_run_id()?)?;

        eprintln!(
            "{} {}",
            ui::muted_stderr("Run artifacts:"),
            self.run_dir.display()
        );
        eprintln!(
            "{} {}",
            ui::muted_stderr("Benchmark database:"),
            self.options.db.display()
        );
        Ok(())
    }

    fn run_tool(&self, subset: &SubsetMetadata, tool: Tool, execution_order: usize) -> Result<()> {
        let artifacts = ToolRunArtifacts::new(&self.run_dir, subset, tool);
        fs::create_dir_all(&artifacts.out_dir)?;
        let command =
            build_tool_command(&self.paths, self.options, subset, tool, &artifacts.out_dir);
        let command_string = artifacts.command_string(&command);
        fs::write(artifacts.command_script(), format!("{command_string}\n"))?;

        eprintln!(
            "{} {} {} {}",
            ui::muted_stderr(format!("[{}]", self.run_id)),
            ui::dataset_stderr(&subset.name),
            ui::muted_stderr(format!("{} pairs", subset.read_pairs)),
            ui::tool_stderr(tool)
        );
        run_command(
            Command::new(&self.paths.hyperfine)
                .arg("--runs")
                .arg(self.options.replicates.to_string())
                .arg("--warmup")
                .arg("1")
                .arg("--export-json")
                .arg(&artifacts.hyperfine_json)
                .arg("--command-name")
                .arg(tool.to_string())
                .arg(&command_string),
        )?;

        let report: HyperfineReport =
            serde_json::from_reader(File::open(&artifacts.hyperfine_json)?)?;
        let result = report
            .results
            .first()
            .ok_or_else(|| color_eyre::eyre::eyre!("hyperfine report had no results"))?;
        let merged_reads = fastq_record_count(&command.merged_output).wrap_err_with(|| {
            format!(
                "failed to count merged FASTQ records for {}: {}",
                tool,
                command.merged_output.display()
            )
        })?;
        validate_tool_run(
            command.tool,
            self.options.mode,
            subset,
            &artifacts.out_dir,
            &command.merged_output,
            merged_reads,
        )?;
        let execution = ToolExecutionRecord {
            run_key: &self.run_key,
            subset,
            tool,
            command_string: &command_string,
            execution_order,
            result,
            merged_reads,
            expected_merged_read_records: merged_reads,
            output_dir: &artifacts.out_dir,
        };
        let artifact_records = artifacts.artifacts(&command);
        self.database.replace_tool_evidence(
            &execution,
            &artifact_records
                .iter()
                .map(|artifact| ArtifactRecord {
                    kind: artifact.kind,
                    path: &artifact.path,
                    requirement: artifact.requirement,
                })
                .collect::<Vec<_>>(),
            &command.merged_output,
        )
    }
}

struct ToolRunArtifacts {
    tool: Tool,
    out_dir: std::path::PathBuf,
    stdout_log: std::path::PathBuf,
    stderr_log: std::path::PathBuf,
    hyperfine_json: std::path::PathBuf,
}

impl ToolRunArtifacts {
    fn new(run_dir: &Path, subset: &SubsetMetadata, tool: Tool) -> Self {
        let out_dir = run_dir
            .join(&subset.name)
            .join(format!("{}_pairs", subset.read_pairs))
            .join(tool.to_string());
        Self {
            stdout_log: out_dir.join(format!("{tool}.stdout.log")),
            stderr_log: out_dir.join(format!("{tool}.stderr.log")),
            hyperfine_json: out_dir.join("hyperfine.json"),
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

    fn artifacts(&self, command: &ToolCommand) -> Vec<ToolArtifact> {
        let mut artifacts = vec![
            ToolArtifact::required(ArtifactKind::Command, self.command_script()),
            ToolArtifact::required(ArtifactKind::StdoutLog, self.stdout_log.clone()),
            ToolArtifact::required(ArtifactKind::StderrLog, self.stderr_log.clone()),
            ToolArtifact::required(ArtifactKind::HyperfineJson, self.hyperfine_json.clone()),
            ToolArtifact::required(ArtifactKind::MergedFastq, command.merged_output.clone()),
        ];

        artifacts.extend(self.tool_artifacts());
        artifacts
    }

    fn tool_artifacts(&self) -> Vec<ToolArtifact> {
        match self.tool {
            Tool::Pairasm => vec![
                ToolArtifact::required(
                    ArtifactKind::PairasmSummaryJson,
                    self.out_dir.join("pairasm.summary.json"),
                ),
                ToolArtifact::optional(
                    ArtifactKind::PairasmUnmergedFastq,
                    self.out_dir.join("pairasm.unmerged.fastq"),
                ),
            ],
            Tool::Fastp => vec![
                ToolArtifact::required(ArtifactKind::FastpJson, self.out_dir.join("fastp.json")),
                ToolArtifact::required(ArtifactKind::FastpHtml, self.out_dir.join("fastp.html")),
                ToolArtifact::optional(
                    ArtifactKind::FastpUnpaired1Fastq,
                    self.out_dir.join("fastp.unpaired1.fastq"),
                ),
                ToolArtifact::optional(
                    ArtifactKind::FastpUnpaired2Fastq,
                    self.out_dir.join("fastp.unpaired2.fastq"),
                ),
                ToolArtifact::optional(
                    ArtifactKind::FastpFailedFastq,
                    self.out_dir.join("fastp.failed.fastq"),
                ),
            ],
            Tool::Bbmerge => vec![
                ToolArtifact::optional(
                    ArtifactKind::BbmergeUnmerged1Fastq,
                    self.out_dir.join("bbmerge.unmerged1.fastq"),
                ),
                ToolArtifact::optional(
                    ArtifactKind::BbmergeUnmerged2Fastq,
                    self.out_dir.join("bbmerge.unmerged2.fastq"),
                ),
            ],
            Tool::Vsearch => vec![
                ToolArtifact::optional(
                    ArtifactKind::VsearchUnmerged1Fastq,
                    self.out_dir.join("vsearch.unmerged1.fastq"),
                ),
                ToolArtifact::optional(
                    ArtifactKind::VsearchUnmerged2Fastq,
                    self.out_dir.join("vsearch.unmerged2.fastq"),
                ),
            ],
        }
    }
}

struct ToolArtifact {
    kind: ArtifactKind,
    path: std::path::PathBuf,
    requirement: ArtifactRequirement,
}

impl ToolArtifact {
    fn required(kind: ArtifactKind, path: std::path::PathBuf) -> Self {
        Self {
            kind,
            path,
            requirement: ArtifactRequirement::Required,
        }
    }

    fn optional(kind: ArtifactKind, path: std::path::PathBuf) -> Self {
        Self {
            kind,
            path,
            requirement: ArtifactRequirement::Optional,
        }
    }
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

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::model::{SubsetMetadata, Tool, ToolCommand};

    use super::ToolRunArtifacts;

    #[test]
    fn tool_run_artifacts_have_stable_kinds_and_paths() {
        let subset = SubsetMetadata {
            name: "dataset-a".to_owned(),
            accession: "DRR000000".to_owned(),
            read_pairs: 100,
            r1: PathBuf::from("r1.fastq.gz"),
            r2: PathBuf::from("r2.fastq.gz"),
        };
        let command = ToolCommand {
            tool: Tool::Pairasm,
            args: Vec::new(),
            merged_output: PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.merged.fastq"),
        };
        let artifacts =
            ToolRunArtifacts::new(PathBuf::from("run").as_path(), &subset, Tool::Pairasm);

        let artifact_specs = artifacts
            .artifacts(&command)
            .into_iter()
            .map(|artifact| (artifact.kind.to_string(), artifact.path))
            .collect::<Vec<_>>();

        assert_eq!(artifact_specs.len(), 7);
        assert!(artifact_specs.contains(&(
            "command".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/command.sh")
        )));
        assert!(artifact_specs.contains(&(
            "stdout_log".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.stdout.log")
        )));
        assert!(artifact_specs.contains(&(
            "stderr_log".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.stderr.log")
        )));
        assert!(artifact_specs.contains(&(
            "hyperfine_json".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/hyperfine.json")
        )));
        assert!(artifact_specs.contains(&(
            "merged_fastq".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.merged.fastq")
        )));
        assert!(artifact_specs.contains(&(
            "pairasm_summary_json".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.summary.json")
        )));
        assert!(artifact_specs.contains(&(
            "pairasm_unmerged_fastq".to_owned(),
            PathBuf::from("run/dataset-a/100_pairs/pairasm/pairasm.unmerged.fastq")
        )));
    }
}
