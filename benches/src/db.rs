use std::{fs, path::Path, process::Command};

use color_eyre::eyre::{Result, WrapErr};
use duckdb::{Connection, params};

use crate::{
    cli::RunOptions,
    model::{HyperfineResult, SubsetMetadata, ToolCommand, ToolPaths},
};

pub struct BenchmarkDb {
    connection: Connection,
}

pub struct RunRecord<'a> {
    pub run_key: &'a str,
    pub run_label: &'a str,
    pub created_at: &'a str,
    pub run_dir: &'a Path,
    pub options: &'a RunOptions,
    pub vcs: &'a VcsMetadata,
}

pub struct ToolExecutionRecord<'a> {
    pub run_key: &'a str,
    pub subset: &'a SubsetMetadata,
    pub command: &'a ToolCommand,
    pub command_string: &'a str,
    pub result: &'a HyperfineResult,
    pub merged_reads: usize,
    pub manifest_rows: usize,
    pub output_dir: &'a Path,
}

#[derive(Debug, Default)]
pub struct VcsMetadata {
    pub vcs_kind: Option<String>,
    pub change_id: Option<String>,
    pub commit_id: Option<String>,
    pub description: Option<String>,
    pub working_copy_dirty: Option<bool>,
}

impl BenchmarkDb {
    pub fn open(path: &Path) -> Result<Self> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let connection = Connection::open(path)
            .wrap_err_with(|| format!("failed to open benchmark database {}", path.display()))?;
        let db = Self { connection };
        db.migrate()?;
        Ok(db)
    }

    pub fn insert_run(&self, record: &RunRecord<'_>) -> Result<()> {
        self.connection.execute(
            "INSERT INTO benchmark_runs (
                run_key, run_label, created_at, benchmark_mode, read_pairs, replicates,
                threads, output_compression, run_dir, config_path, data_root, vcs_kind,
                vcs_change_id, vcs_commit_id, vcs_description, vcs_working_copy_dirty
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16)",
            params![
                record.run_key,
                record.run_label,
                record.created_at,
                record.options.mode.to_string(),
                record.options.read_pairs as i64,
                record.options.replicates as i64,
                record.options.threads as i64,
                record.options.output_compression.to_string(),
                record.run_dir.to_string_lossy().to_string(),
                record.options.common.config.to_string_lossy().to_string(),
                record
                    .options
                    .common
                    .data_root
                    .to_string_lossy()
                    .to_string(),
                record.vcs.vcs_kind.as_deref(),
                record.vcs.change_id.as_deref(),
                record.vcs.commit_id.as_deref(),
                record.vcs.description.as_deref(),
                record.vcs.working_copy_dirty.map(i64::from),
            ],
        )?;
        Ok(())
    }

    pub fn insert_dataset_subset(&self, run_key: &str, subset: &SubsetMetadata) -> Result<()> {
        self.connection.execute(
            "INSERT INTO dataset_subsets (
                run_key, dataset_name, accession, read_pairs, r1_path, r2_path, r1_bytes, r2_bytes
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                run_key,
                subset.name,
                subset.accession,
                subset.read_pairs as i64,
                subset.r1.to_string_lossy().to_string(),
                subset.r2.to_string_lossy().to_string(),
                file_size(&subset.r1)? as i64,
                file_size(&subset.r2)? as i64,
            ],
        )?;
        Ok(())
    }

    pub fn insert_tool_versions(&self, run_key: &str, paths: &ToolPaths) -> Result<()> {
        for (tool, path) in [
            ("pairasm", &paths.pairasm),
            ("fastp", &paths.fastp),
            ("bbmerge", &paths.bbmerge),
            ("vsearch", &paths.vsearch),
            ("hyperfine", &paths.hyperfine),
        ] {
            self.connection.execute(
                "INSERT INTO tool_versions (run_key, tool, path) VALUES (?1, ?2, ?3)",
                params![run_key, tool, path.to_string_lossy().to_string()],
            )?;
        }
        Ok(())
    }

    pub fn insert_tool_execution(&self, record: &ToolExecutionRecord<'_>) -> Result<()> {
        self.connection.execute(
            "INSERT INTO tool_executions (
                run_key, dataset_name, tool, command, hyperfine_mean_s, hyperfine_median_s,
                hyperfine_stddev_s, hyperfine_min_s, hyperfine_max_s, hyperfine_user_s,
                hyperfine_system_s, merged_reads, manifest_rows, output_dir
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14)",
            params![
                record.run_key,
                record.subset.name,
                record.command.tool.name(),
                record.command_string,
                record.result.mean,
                record.result.median,
                record.result.stddev,
                record.result.min,
                record.result.max,
                record.result.user,
                record.result.system,
                record.merged_reads as i64,
                record.manifest_rows as i64,
                record.output_dir.to_string_lossy().to_string(),
            ],
        )?;
        Ok(())
    }

    pub fn insert_artifact(
        &self,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
        kind: &str,
        path: &Path,
    ) -> Result<()> {
        let bytes = path.exists().then(|| file_size(path)).transpose()?;
        self.connection.execute(
            "INSERT INTO artifacts (run_key, dataset_name, tool, artifact_kind, path, bytes)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![
                run_key,
                dataset_name,
                tool,
                kind,
                path.to_string_lossy().to_string(),
                bytes.map(|value| value as i64),
            ],
        )?;
        Ok(())
    }

    fn migrate(&self) -> Result<()> {
        self.connection.execute_batch(SCHEMA)?;
        Ok(())
    }
}

pub fn collect_vcs_metadata() -> VcsMetadata {
    collect_jj_metadata()
        .or_else(collect_git_metadata)
        .unwrap_or_default()
}

fn collect_jj_metadata() -> Option<VcsMetadata> {
    let status = command_stdout(Command::new("jj").arg("status"))?;
    let log = command_stdout(Command::new("jj").arg("log").args([
        "-r",
        "@",
        "--no-graph",
        "--limit",
        "1",
    ]))?;
    Some(VcsMetadata {
        vcs_kind: Some("jj".to_owned()),
        change_id: parse_jj_change_id(&log),
        commit_id: parse_jj_commit_id(&log),
        description: parse_jj_description(&log),
        working_copy_dirty: Some(status.contains("Working copy changes:")),
    })
}

fn collect_git_metadata() -> Option<VcsMetadata> {
    let commit_id = command_stdout(Command::new("git").args(["rev-parse", "HEAD"]))?;
    let description = command_stdout(Command::new("git").args(["log", "-1", "--format=%s"]))?;
    let status = command_stdout(Command::new("git").args(["status", "--porcelain"]))?;
    Some(VcsMetadata {
        vcs_kind: Some("git".to_owned()),
        change_id: None,
        commit_id: Some(commit_id.trim().to_owned()),
        description: Some(description.trim().to_owned()),
        working_copy_dirty: Some(!status.trim().is_empty()),
    })
}

fn command_stdout(command: &mut Command) -> Option<String> {
    let output = command.output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8(output.stdout).ok())?
}

fn parse_jj_change_id(log: &str) -> Option<String> {
    log.lines()
        .find(|line| line.contains(' '))
        .and_then(|line| line.split_whitespace().next())
        .map(str::to_owned)
}

fn parse_jj_commit_id(log: &str) -> Option<String> {
    log.lines()
        .find(|line| line.contains(' '))
        .and_then(|line| line.split_whitespace().nth(3))
        .map(str::to_owned)
}

fn parse_jj_description(log: &str) -> Option<String> {
    log.lines()
        .map(str::trim)
        .find(|line| !line.is_empty() && !line.contains('@'))
        .map(str::to_owned)
}

fn file_size(path: &Path) -> Result<u64> {
    Ok(path.metadata()?.len())
}

const SCHEMA: &str = r"
CREATE TABLE IF NOT EXISTS benchmark_runs (
    run_key TEXT PRIMARY KEY,
    run_label TEXT NOT NULL,
    created_at TEXT NOT NULL,
    benchmark_mode TEXT NOT NULL,
    read_pairs BIGINT NOT NULL,
    replicates BIGINT NOT NULL,
    threads BIGINT NOT NULL,
    output_compression TEXT NOT NULL,
    run_dir TEXT NOT NULL,
    config_path TEXT NOT NULL,
    data_root TEXT NOT NULL,
    vcs_kind TEXT,
    vcs_change_id TEXT,
    vcs_commit_id TEXT,
    vcs_description TEXT,
    vcs_working_copy_dirty BOOLEAN
);

CREATE TABLE IF NOT EXISTS dataset_subsets (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    accession TEXT NOT NULL,
    read_pairs BIGINT NOT NULL,
    r1_path TEXT NOT NULL,
    r2_path TEXT NOT NULL,
    r1_bytes BIGINT NOT NULL,
    r2_bytes BIGINT NOT NULL,
    PRIMARY KEY (run_key, dataset_name)
);

CREATE TABLE IF NOT EXISTS tool_versions (
    run_key TEXT NOT NULL,
    tool TEXT NOT NULL,
    path TEXT NOT NULL,
    version TEXT,
    PRIMARY KEY (run_key, tool)
);

CREATE TABLE IF NOT EXISTS tool_executions (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    tool TEXT NOT NULL,
    command TEXT NOT NULL,
    hyperfine_mean_s DOUBLE NOT NULL,
    hyperfine_median_s DOUBLE NOT NULL,
    hyperfine_stddev_s DOUBLE,
    hyperfine_min_s DOUBLE NOT NULL,
    hyperfine_max_s DOUBLE NOT NULL,
    hyperfine_user_s DOUBLE NOT NULL,
    hyperfine_system_s DOUBLE NOT NULL,
    merged_reads BIGINT NOT NULL,
    manifest_rows BIGINT NOT NULL,
    output_dir TEXT NOT NULL,
    PRIMARY KEY (run_key, dataset_name, tool)
);

CREATE TABLE IF NOT EXISTS artifacts (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    tool TEXT NOT NULL,
    artifact_kind TEXT NOT NULL,
    path TEXT NOT NULL,
    bytes BIGINT,
    PRIMARY KEY (run_key, dataset_name, tool, artifact_kind)
);
";
