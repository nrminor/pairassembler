use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use color_eyre::eyre::{Result, WrapErr};
use duckdb::{Connection, params};

use crate::{
    cli::RunOptions,
    model::{HyperfineResult, SubsetMetadata, ToolCommand, ToolPaths},
    products::MergedProduct,
};

pub struct BenchmarkDb {
    connection: Connection,
    path: PathBuf,
}

pub struct AgreementRow {
    pub run_label: String,
    pub dataset_name: String,
    pub left_tool: String,
    pub right_tool: String,
    pub left_merged: i64,
    pub right_merged: i64,
    pub shared_merged: i64,
    pub left_only: i64,
    pub right_only: i64,
    pub jaccard: f64,
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
    pub merged_product_rows: usize,
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
        let db = Self {
            connection,
            path: path.to_owned(),
        };
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
                record.merged_product_rows as i64,
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

    pub fn insert_merged_products(
        &self,
        run_key: &str,
        subset: &SubsetMetadata,
        command: &ToolCommand,
        products: &[MergedProduct],
    ) -> Result<()> {
        self.connection.execute_batch("BEGIN TRANSACTION")?;
        let insert_result = (|| -> Result<()> {
            let mut statement = self.connection.prepare(
                "INSERT INTO merged_products (
                    run_key, dataset_name, tool, read_id, output_header, merged_len,
                    avg_qual, min_qual, max_qual, sequence_hash, quality_hash
                ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            )?;
            for product in products {
                statement.execute(params![
                    run_key,
                    subset.name,
                    command.tool.name(),
                    product.read_id.as_str(),
                    product.output_header.as_str(),
                    product.merged_len as i64,
                    product.avg_qual,
                    i64::from(product.min_qual),
                    i64::from(product.max_qual),
                    product.sequence_hash.as_str(),
                    product.quality_hash.as_str(),
                ])?;
            }
            Ok(())
        })();

        if let Err(error) = insert_result {
            self.connection.execute_batch("ROLLBACK")?;
            return Err(error);
        }

        self.connection.execute_batch("COMMIT")?;
        Ok(())
    }

    pub fn latest_run_label(&self) -> Result<String> {
        let mut statement = self.connection.prepare(
            "SELECT run_label FROM benchmark_runs ORDER BY created_at DESC, run_label DESC LIMIT 1",
        )?;
        let mut rows = statement.query([])?;
        if let Some(row) = rows.next()? {
            return Ok(row.get(0)?);
        }

        color_eyre::eyre::bail!(
            "no benchmark runs have been recorded yet\n\nTo run the standard real-data comparison:\n\n  just benchmark\n\nThis checks external tools, fetches configured ENA inputs, prepares deterministic subsets, runs each merge tool, validates the outputs, and prints a report.\n\nFor tuned/comparability mode:\n\n  just benchmark-tuned\n\nResults store: {}",
            self.path.display()
        )
    }

    pub fn agreement_rows(&self, run_label: &str) -> Result<Vec<AgreementRow>> {
        let mut statement = self.connection.prepare(
            "WITH selected_run AS (
                SELECT run_key, run_label
                FROM benchmark_runs
                WHERE run_label = ?1
            ),
            product_tools AS (
                SELECT DISTINCT products.run_key, products.dataset_name, products.tool
                FROM merged_products products
                JOIN selected_run ON selected_run.run_key = products.run_key
            ),
            pairs AS (
                SELECT
                    left_tool.run_key,
                    left_tool.dataset_name,
                    left_tool.tool AS left_tool,
                    right_tool.tool AS right_tool
                FROM product_tools left_tool
                JOIN product_tools right_tool
                    ON right_tool.run_key = left_tool.run_key
                    AND right_tool.dataset_name = left_tool.dataset_name
                    AND right_tool.tool > left_tool.tool
            ),
            tool_counts AS (
                SELECT products.run_key, products.dataset_name, products.tool, COUNT(*) AS merged
                FROM merged_products products
                JOIN selected_run ON selected_run.run_key = products.run_key
                GROUP BY products.run_key, products.dataset_name, products.tool
            ),
            intersections AS (
                SELECT
                    left_product.run_key,
                    left_product.dataset_name,
                    left_product.tool AS left_tool,
                    right_product.tool AS right_tool,
                    COUNT(*) AS shared_merged
                FROM merged_products left_product
                JOIN merged_products right_product
                    ON right_product.run_key = left_product.run_key
                    AND right_product.dataset_name = left_product.dataset_name
                    AND right_product.read_id = left_product.read_id
                    AND right_product.tool > left_product.tool
                JOIN selected_run ON selected_run.run_key = left_product.run_key
                GROUP BY
                    left_product.run_key,
                    left_product.dataset_name,
                    left_product.tool,
                    right_product.tool
            )
            SELECT
                selected_run.run_label,
                pairs.dataset_name,
                pairs.left_tool,
                pairs.right_tool,
                left_counts.merged AS left_merged,
                right_counts.merged AS right_merged,
                COALESCE(intersections.shared_merged, 0) AS shared_merged,
                left_counts.merged - COALESCE(intersections.shared_merged, 0) AS left_only,
                right_counts.merged - COALESCE(intersections.shared_merged, 0) AS right_only,
                round(
                    CASE
                        WHEN left_counts.merged + right_counts.merged - COALESCE(intersections.shared_merged, 0) = 0
                        THEN 0.0
                        ELSE
                            COALESCE(intersections.shared_merged, 0)::DOUBLE
                            / (left_counts.merged + right_counts.merged - COALESCE(intersections.shared_merged, 0))::DOUBLE
                    END,
                    6
                ) AS jaccard
            FROM pairs
            JOIN selected_run ON selected_run.run_key = pairs.run_key
            JOIN tool_counts left_counts
                ON left_counts.run_key = pairs.run_key
                AND left_counts.dataset_name = pairs.dataset_name
                AND left_counts.tool = pairs.left_tool
            JOIN tool_counts right_counts
                ON right_counts.run_key = pairs.run_key
                AND right_counts.dataset_name = pairs.dataset_name
                AND right_counts.tool = pairs.right_tool
            LEFT JOIN intersections
                ON intersections.run_key = pairs.run_key
                AND intersections.dataset_name = pairs.dataset_name
                AND intersections.left_tool = pairs.left_tool
                AND intersections.right_tool = pairs.right_tool
            ORDER BY pairs.dataset_name, pairs.left_tool, pairs.right_tool",
        )?;

        let rows = statement.query_map(params![run_label], |row| {
            Ok(AgreementRow {
                run_label: row.get(0)?,
                dataset_name: row.get(1)?,
                left_tool: row.get(2)?,
                right_tool: row.get(3)?,
                left_merged: row.get(4)?,
                right_merged: row.get(5)?,
                shared_merged: row.get(6)?,
                left_only: row.get(7)?,
                right_only: row.get(8)?,
                jaccard: row.get(9)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
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

CREATE TABLE IF NOT EXISTS merged_products (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    tool TEXT NOT NULL,
    read_id TEXT NOT NULL,
    output_header TEXT NOT NULL,
    merged_len BIGINT NOT NULL,
    avg_qual DOUBLE NOT NULL,
    min_qual BIGINT NOT NULL,
    max_qual BIGINT NOT NULL,
    sequence_hash TEXT NOT NULL,
    quality_hash TEXT NOT NULL,
    PRIMARY KEY (run_key, dataset_name, tool, read_id)
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
