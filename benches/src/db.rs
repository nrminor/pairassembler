use std::{
    fs,
    path::{Path, PathBuf},
};

use color_eyre::eyre::{Result, WrapErr, bail};
use duckdb::{Connection, params};

use crate::{
    artifacts::ArtifactRecord,
    cli::RunOptions,
    config::version_string,
    config::{SubsetMetadata, ToolPaths},
    products::for_each_merged_read_record,
    tool::Tool,
    vcs::VcsMetadata,
};

pub(crate) struct BenchmarkDb {
    connection: Connection,
    path: PathBuf,
}

pub(crate) struct AgreementRow {
    pub(crate) dataset_name: String,
    pub(crate) left_tool: String,
    pub(crate) right_tool: String,
    pub(crate) left_merged: i64,
    pub(crate) right_merged: i64,
    pub(crate) shared_merged: i64,
    pub(crate) left_only: i64,
    pub(crate) right_only: i64,
    pub(crate) jaccard: f64,
}

pub(crate) struct RunSummary {
    pub(crate) run_key: String,
    pub(crate) benchmark_mode: String,
    pub(crate) read_pairs: i64,
    pub(crate) replicates: i64,
    pub(crate) threads: i64,
    pub(crate) run_dir: String,
    pub(crate) status: String,
    pub(crate) completed_at: Option<String>,
}

pub(crate) struct ToolResultRow {
    pub(crate) run_key: String,
    pub(crate) benchmark_mode: String,
    pub(crate) dataset_name: String,
    pub(crate) accession: String,
    pub(crate) read_pairs: i64,
    pub(crate) tool: String,
    pub(crate) replicates: i64,
    pub(crate) threads: i64,
    pub(crate) mean_s: f64,
    pub(crate) median_s: f64,
    pub(crate) stddev_s: Option<f64>,
    pub(crate) min_s: f64,
    pub(crate) max_s: f64,
    pub(crate) user_s: f64,
    pub(crate) system_s: f64,
    pub(crate) merged_reads: i64,
    pub(crate) merged_read_records: i64,
    pub(crate) r1_bytes: i64,
    pub(crate) r2_bytes: i64,
    pub(crate) output_dir: String,
}

pub(crate) struct RunRecord<'a> {
    pub(crate) run_key: &'a str,
    pub(crate) created_at: &'a str,
    pub(crate) run_dir: &'a Path,
    pub(crate) options: &'a RunOptions,
    pub(crate) vcs: &'a VcsMetadata,
}

pub(crate) struct ToolExecutionRecord<'a> {
    pub(crate) run_key: &'a str,
    pub(crate) subset: &'a SubsetMetadata,
    pub(crate) tool: Tool,
    pub(crate) command_string: &'a str,
    pub(crate) execution_order: usize,
    pub(crate) timing: &'a TimingRecord,
    pub(crate) merged_reads: usize,
    pub(crate) expected_merged_read_records: usize,
    pub(crate) output_dir: &'a Path,
}

pub(crate) struct TimingRecord {
    pub(crate) mean_s: f64,
    pub(crate) median_s: f64,
    pub(crate) stddev_s: Option<f64>,
    pub(crate) min_s: f64,
    pub(crate) max_s: f64,
    pub(crate) user_s: f64,
    pub(crate) system_s: f64,
}

impl BenchmarkDb {
    pub(crate) fn open(path: &Path) -> Result<Self> {
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

    pub(crate) fn open_existing(path: &Path) -> Result<Self> {
        if !path.exists() {
            color_eyre::eyre::bail!(
                "benchmark database does not exist: {}\n\nRun a benchmark first or pass the correct --db path.",
                path.display()
            );
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

    pub(crate) fn insert_run(&self, record: &RunRecord<'_>) -> Result<()> {
        self.connection.execute(
            "INSERT INTO benchmark_runs (
                run_key, created_at, benchmark_mode, read_pairs, replicates,
                threads, run_dir, config_path, data_root, vcs_kind,
                vcs_change_id, vcs_commit_id, vcs_description, vcs_working_copy_dirty
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14)",
            params![
                record.run_key,
                record.created_at,
                record.options.mode.to_string(),
                record.options.read_pairs as i64,
                record.options.replicates as i64,
                record.options.threads as i64,
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

    pub(crate) fn mark_run_completed(&self, run_key: &str, completed_at: &str) -> Result<()> {
        self.connection.execute(
            "UPDATE benchmark_runs SET status = 'completed', completed_at = ?2 WHERE run_key = ?1",
            params![run_key, completed_at],
        )?;
        Ok(())
    }

    pub(crate) fn insert_dataset_subset(
        &self,
        run_key: &str,
        subset: &SubsetMetadata,
    ) -> Result<()> {
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

    pub(crate) fn insert_tool_versions(&self, run_key: &str, paths: &ToolPaths) -> Result<()> {
        for (tool, path) in paths.versioned_tools() {
            let version = version_string(path).ok();
            self.connection.execute(
                "INSERT INTO tool_versions (run_key, tool, path, version) VALUES (?1, ?2, ?3, ?4)",
                params![run_key, tool, path.to_string_lossy().to_string(), version],
            )?;
        }
        Ok(())
    }

    fn record_tool_execution(&self, record: &ToolExecutionRecord<'_>) -> Result<()> {
        let tool = record.tool.to_string();
        let output_dir = record.output_dir.to_string_lossy().to_string();
        let updated = self.connection.execute(
            "UPDATE tool_executions SET
                command = ?4,
                execution_order = ?5,
                hyperfine_mean_s = ?6,
                hyperfine_median_s = ?7,
                hyperfine_stddev_s = ?8,
                hyperfine_min_s = ?9,
                hyperfine_max_s = ?10,
                hyperfine_user_s = ?11,
                hyperfine_system_s = ?12,
                merged_reads = ?13,
                merged_read_records = ?14,
                output_dir = ?15
             WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
            params![
                record.run_key,
                record.subset.name,
                tool,
                record.command_string,
                record.execution_order as i64,
                record.timing.mean_s,
                record.timing.median_s,
                record.timing.stddev_s,
                record.timing.min_s,
                record.timing.max_s,
                record.timing.user_s,
                record.timing.system_s,
                record.merged_reads as i64,
                record.expected_merged_read_records as i64,
                output_dir,
            ],
        )?;
        if updated > 0 {
            return Ok(());
        }

        self.connection.execute(
            "INSERT INTO tool_executions (
                run_key, dataset_name, tool, command, execution_order, hyperfine_mean_s, hyperfine_median_s,
                hyperfine_stddev_s, hyperfine_min_s, hyperfine_max_s, hyperfine_user_s,
                hyperfine_system_s, merged_reads, merged_read_records, output_dir
            ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15)",
            params![
                record.run_key,
                record.subset.name,
                tool,
                record.command_string,
                record.execution_order as i64,
                record.timing.mean_s,
                record.timing.median_s,
                record.timing.stddev_s,
                record.timing.min_s,
                record.timing.max_s,
                record.timing.user_s,
                record.timing.system_s,
                record.merged_reads as i64,
                record.expected_merged_read_records as i64,
                output_dir,
            ],
        )?;
        Ok(())
    }

    fn insert_artifact(
        &self,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
        artifact: &ArtifactRecord<'_>,
    ) -> Result<()> {
        let bytes = artifact
            .path
            .exists()
            .then(|| file_size(artifact.path))
            .transpose()?;
        self.connection.execute(
            "INSERT INTO artifacts (run_key, dataset_name, tool, artifact_kind, path, required, bytes)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
            params![
                run_key,
                dataset_name,
                tool,
                artifact.kind.to_string(),
                artifact.path.to_string_lossy().to_string(),
                i64::from(artifact.requirement.is_required()),
                bytes.map(|value| value as i64),
            ],
        )?;
        Ok(())
    }

    pub(crate) fn replace_tool_evidence(
        &self,
        execution: &ToolExecutionRecord<'_>,
        artifacts: &[ArtifactRecord<'_>],
        merged_fastq: &Path,
    ) -> Result<()> {
        self.connection.execute_batch("BEGIN TRANSACTION")?;
        let result = (|| -> Result<()> {
            let tool = execution.tool.to_string();
            self.connection.execute(
                "DELETE FROM artifacts WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
                params![execution.run_key, execution.subset.name, tool.as_str()],
            )?;
            self.connection.execute(
                "DELETE FROM merged_read_records WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
                params![execution.run_key, execution.subset.name, tool.as_str()],
            )?;
            self.record_tool_execution(execution)?;
            for artifact in artifacts {
                self.insert_artifact(
                    execution.run_key,
                    &execution.subset.name,
                    tool.as_str(),
                    artifact,
                )?;
            }
            let merged_read_records = self.insert_merged_read_records_from_fastq(
                execution.run_key,
                execution.subset,
                execution.tool,
                merged_fastq,
            )?;
            if merged_read_records != execution.expected_merged_read_records {
                bail!(
                    "merged read count changed while recording {} {} {}: counted {}, inserted {}",
                    execution.subset.name,
                    execution.tool,
                    merged_fastq.display(),
                    execution.expected_merged_read_records,
                    merged_read_records
                );
            }
            Ok(())
        })();

        if let Err(error) = result {
            self.connection.execute_batch("ROLLBACK")?;
            return Err(error);
        }

        self.connection.execute_batch("COMMIT")?;
        Ok(())
    }

    fn insert_merged_read_records_from_fastq(
        &self,
        run_key: &str,
        subset: &SubsetMetadata,
        tool: Tool,
        merged_fastq: &Path,
    ) -> Result<usize> {
        const MERGED_READ_RECORD_COLUMNS: &[&str] = &[
            "run_key",
            "dataset_name",
            "tool",
            "read_id",
            "output_header",
            "merged_len",
            "avg_qual",
            "min_qual",
            "max_qual",
            "sequence_hash",
            "quality_hash",
        ];

        let tool = tool.to_string();
        let dataset_name = subset.name.as_str();
        let count = {
            let mut appender = self
                .connection
                .appender_with_columns("merged_read_records", MERGED_READ_RECORD_COLUMNS)?;
            let count = for_each_merged_read_record(merged_fastq, |merged_read| {
                appender.append_row(params![
                    run_key,
                    dataset_name,
                    tool.as_str(),
                    merged_read.read_id.as_str(),
                    merged_read.output_header.as_str(),
                    merged_read.merged_len as i64,
                    merged_read.avg_qual,
                    i64::from(merged_read.min_qual),
                    i64::from(merged_read.max_qual),
                    merged_read.sequence_hash.as_str(),
                    merged_read.quality_hash.as_str(),
                ])?;
                Ok(())
            })?;
            appender.flush()?;
            count
        };

        Ok(count)
    }

    pub(crate) fn latest_completed_run_key_for_mode(&self, benchmark_mode: &str) -> Result<String> {
        let mut statement = self.connection.prepare(
            "SELECT run_key
             FROM benchmark_runs
             WHERE benchmark_mode = ?1
             AND status = 'completed'
             ORDER BY created_at DESC, run_key DESC
             LIMIT 1",
        )?;
        let mut rows = statement.query(params![benchmark_mode])?;
        if let Some(row) = rows.next()? {
            return Ok(row.get(0)?);
        }

        color_eyre::eyre::bail!(
            "no {benchmark_mode} benchmark runs have been recorded yet\n\nTo run the standard real-data comparison:\n\n  just benchmark\n\nThis checks external tools, fetches configured ENA inputs, prepares deterministic subsets, runs each merge tool, validates the outputs, and prints a report.\n\nFor tuned/comparability mode:\n\n  just benchmark-tuned\n\nResults store: {}",
            self.path.display()
        )
    }

    pub(crate) fn run_summary(&self, run_key: &str) -> Result<RunSummary> {
        let mut statement = self.connection.prepare(
            "SELECT run_key, benchmark_mode, read_pairs, replicates, threads, run_dir, status, completed_at
             FROM benchmark_runs
             WHERE run_key = ?1
             ORDER BY created_at DESC
             LIMIT 1",
        )?;
        let mut rows = statement.query(params![run_key])?;
        if let Some(row) = rows.next()? {
            return Ok(RunSummary {
                run_key: row.get(0)?,
                benchmark_mode: row.get(1)?,
                read_pairs: row.get(2)?,
                replicates: row.get(3)?,
                threads: row.get(4)?,
                run_dir: row.get(5)?,
                status: row.get(6)?,
                completed_at: row.get(7)?,
            });
        }

        color_eyre::eyre::bail!("no benchmark run found for key {run_key}")
    }

    pub(crate) fn agreement_rows(&self, run_key: &str) -> Result<Vec<AgreementRow>> {
        let mut statement = self.connection.prepare(
            "WITH selected_run AS (
                SELECT run_key
                FROM benchmark_runs
                WHERE run_key = ?1
                AND status = 'completed'
            ),
            execution_tools AS (
                SELECT DISTINCT executions.run_key, executions.dataset_name, executions.tool
                FROM tool_executions executions
                JOIN selected_run ON selected_run.run_key = executions.run_key
            ),
            pairs AS (
                SELECT
                    left_tool.run_key,
                    left_tool.dataset_name,
                    left_tool.tool AS left_tool,
                    right_tool.tool AS right_tool
                FROM execution_tools left_tool
                JOIN execution_tools right_tool
                    ON right_tool.run_key = left_tool.run_key
                    AND right_tool.dataset_name = left_tool.dataset_name
                    AND right_tool.tool > left_tool.tool
            ),
            tool_counts AS (
                SELECT records.run_key, records.dataset_name, records.tool, COUNT(*) AS merged
                FROM merged_read_records records
                JOIN selected_run ON selected_run.run_key = records.run_key
                GROUP BY records.run_key, records.dataset_name, records.tool
            ),
            intersections AS (
                SELECT
                    left_record.run_key,
                    left_record.dataset_name,
                    left_record.tool AS left_tool,
                    right_record.tool AS right_tool,
                    COUNT(*) AS shared_merged
                FROM merged_read_records left_record
                JOIN merged_read_records right_record
                    ON right_record.run_key = left_record.run_key
                    AND right_record.dataset_name = left_record.dataset_name
                    AND right_record.read_id = left_record.read_id
                    AND right_record.tool > left_record.tool
                JOIN selected_run ON selected_run.run_key = left_record.run_key
                GROUP BY
                    left_record.run_key,
                    left_record.dataset_name,
                    left_record.tool,
                    right_record.tool
            )
            SELECT
                pairs.dataset_name,
                pairs.left_tool,
                pairs.right_tool,
                COALESCE(left_counts.merged, 0) AS left_merged,
                COALESCE(right_counts.merged, 0) AS right_merged,
                COALESCE(intersections.shared_merged, 0) AS shared_merged,
                COALESCE(left_counts.merged, 0) - COALESCE(intersections.shared_merged, 0) AS left_only,
                COALESCE(right_counts.merged, 0) - COALESCE(intersections.shared_merged, 0) AS right_only,
                round(
                    CASE
                        WHEN COALESCE(left_counts.merged, 0) + COALESCE(right_counts.merged, 0) - COALESCE(intersections.shared_merged, 0) = 0
                        THEN 0.0
                        ELSE
                            COALESCE(intersections.shared_merged, 0)::DOUBLE
                            / (COALESCE(left_counts.merged, 0) + COALESCE(right_counts.merged, 0) - COALESCE(intersections.shared_merged, 0))::DOUBLE
                    END,
                    6
                ) AS jaccard
            FROM pairs
            JOIN selected_run ON selected_run.run_key = pairs.run_key
            LEFT JOIN tool_counts left_counts
                ON left_counts.run_key = pairs.run_key
                AND left_counts.dataset_name = pairs.dataset_name
                AND left_counts.tool = pairs.left_tool
            LEFT JOIN tool_counts right_counts
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

        let rows = statement.query_map(params![run_key], |row| {
            Ok(AgreementRow {
                dataset_name: row.get(0)?,
                left_tool: row.get(1)?,
                right_tool: row.get(2)?,
                left_merged: row.get(3)?,
                right_merged: row.get(4)?,
                shared_merged: row.get(5)?,
                left_only: row.get(6)?,
                right_only: row.get(7)?,
                jaccard: row.get(8)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    pub(crate) fn tool_result_rows(&self, run_key: &str) -> Result<Vec<ToolResultRow>> {
        let mut statement = self.connection.prepare(
            "SELECT
                runs.run_key,
                runs.benchmark_mode,
                subsets.dataset_name,
                subsets.accession,
                subsets.read_pairs,
                executions.tool,
                runs.replicates,
                runs.threads,
                executions.hyperfine_mean_s,
                executions.hyperfine_median_s,
                executions.hyperfine_stddev_s,
                executions.hyperfine_min_s,
                executions.hyperfine_max_s,
                executions.hyperfine_user_s,
                executions.hyperfine_system_s,
                executions.merged_reads,
                executions.merged_read_records,
                subsets.r1_bytes,
                subsets.r2_bytes,
                executions.output_dir
             FROM benchmark_runs runs
             JOIN dataset_subsets subsets ON subsets.run_key = runs.run_key
             JOIN tool_executions executions
                ON executions.run_key = runs.run_key
                AND executions.dataset_name = subsets.dataset_name
             WHERE runs.run_key = ?1
             AND runs.status = 'completed'
             ORDER BY subsets.dataset_name, executions.tool",
        )?;

        let rows = statement.query_map(params![run_key], |row| {
            Ok(ToolResultRow {
                run_key: row.get(0)?,
                benchmark_mode: row.get(1)?,
                dataset_name: row.get(2)?,
                accession: row.get(3)?,
                read_pairs: row.get(4)?,
                tool: row.get(5)?,
                replicates: row.get(6)?,
                threads: row.get(7)?,
                mean_s: row.get(8)?,
                median_s: row.get(9)?,
                stddev_s: row.get(10)?,
                min_s: row.get(11)?,
                max_s: row.get(12)?,
                user_s: row.get(13)?,
                system_s: row.get(14)?,
                merged_reads: row.get(15)?,
                merged_read_records: row.get(16)?,
                r1_bytes: row.get(17)?,
                r2_bytes: row.get(18)?,
                output_dir: row.get(19)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    fn migrate(&self) -> Result<()> {
        self.connection.execute_batch(SCHEMA)?;
        Ok(())
    }
}

fn file_size(path: &Path) -> Result<u64> {
    Ok(path.metadata()?.len())
}

const SCHEMA: &str = r"
CREATE TABLE IF NOT EXISTS benchmark_runs (
    run_key TEXT PRIMARY KEY,
    created_at TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'running' CHECK (status IN ('running', 'completed')),
    completed_at TEXT,
    benchmark_mode TEXT NOT NULL,
    read_pairs BIGINT NOT NULL,
    replicates BIGINT NOT NULL,
    threads BIGINT NOT NULL,
    run_dir TEXT NOT NULL,
    config_path TEXT NOT NULL,
    data_root TEXT NOT NULL,
    vcs_kind TEXT,
    vcs_change_id TEXT,
    vcs_commit_id TEXT,
    vcs_description TEXT,
    vcs_working_copy_dirty BOOLEAN,
    CHECK ((status = 'completed') = (completed_at IS NOT NULL))
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
    PRIMARY KEY (run_key, dataset_name),
    FOREIGN KEY (run_key) REFERENCES benchmark_runs(run_key)
);

CREATE TABLE IF NOT EXISTS tool_versions (
    run_key TEXT NOT NULL,
    tool TEXT NOT NULL,
    path TEXT NOT NULL,
    version TEXT,
    PRIMARY KEY (run_key, tool),
    FOREIGN KEY (run_key) REFERENCES benchmark_runs(run_key)
);

CREATE TABLE IF NOT EXISTS tool_executions (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    tool TEXT NOT NULL,
    command TEXT NOT NULL,
    execution_order BIGINT NOT NULL,
    hyperfine_mean_s DOUBLE NOT NULL,
    hyperfine_median_s DOUBLE NOT NULL,
    hyperfine_stddev_s DOUBLE,
    hyperfine_min_s DOUBLE NOT NULL,
    hyperfine_max_s DOUBLE NOT NULL,
    hyperfine_user_s DOUBLE NOT NULL,
    hyperfine_system_s DOUBLE NOT NULL,
    merged_reads BIGINT NOT NULL,
    merged_read_records BIGINT NOT NULL,
    output_dir TEXT NOT NULL,
    PRIMARY KEY (run_key, dataset_name, tool),
    FOREIGN KEY (run_key, dataset_name) REFERENCES dataset_subsets(run_key, dataset_name)
);

CREATE TABLE IF NOT EXISTS merged_read_records (
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
    PRIMARY KEY (run_key, dataset_name, tool, read_id),
    FOREIGN KEY (run_key, dataset_name, tool) REFERENCES tool_executions(run_key, dataset_name, tool)
);

CREATE TABLE IF NOT EXISTS artifacts (
    run_key TEXT NOT NULL,
    dataset_name TEXT NOT NULL,
    tool TEXT NOT NULL,
    artifact_kind TEXT NOT NULL,
    path TEXT NOT NULL,
    required BOOLEAN NOT NULL,
    bytes BIGINT,
    PRIMARY KEY (run_key, dataset_name, tool, artifact_kind),
    FOREIGN KEY (run_key, dataset_name, tool) REFERENCES tool_executions(run_key, dataset_name, tool),
    CHECK (artifact_kind IN (
        'command',
        'stdout_log',
        'stderr_log',
        'hyperfine_json',
        'merged_fastq',
        'pairasm_summary_json',
        'pairasm_unmerged_fastq',
        'fastp_json',
        'fastp_html',
        'fastp_unpaired1_fastq',
        'fastp_unpaired2_fastq',
        'fastp_failed_fastq',
        'bbmerge_unmerged1_fastq',
        'bbmerge_unmerged2_fastq',
        'vsearch_unmerged1_fastq',
        'vsearch_unmerged2_fastq'
    ))
);
";

#[cfg(test)]
mod tests {
    use crate::{
        artifacts::{ArtifactKind, ArtifactRecord, ArtifactRequirement},
        tool::Tool,
    };

    use super::*;

    #[test]
    fn agreement_rows_include_tools_with_zero_merged_reads() -> Result<()> {
        let path = temp_db_path("agreement-zero-products");
        let db = BenchmarkDb::open(&path)?;
        insert_run(&db, "run-1", "2026-05-08T00:00:00Z", "completed")?;
        insert_dataset_subset(&db, "run-1", "dataset-a")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "pairasm")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "fastp")?;
        insert_merged_read_record(&db, "run-1", "dataset-a", "pairasm", "read-1")?;

        let rows = db.agreement_rows("run-1")?;

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.dataset_name, "dataset-a");
        assert_eq!(row.left_tool, "fastp");
        assert_eq!(row.right_tool, "pairasm");
        assert_eq!(row.left_merged, 0);
        assert_eq!(row.right_merged, 1);
        assert_eq!(row.shared_merged, 0);
        assert_eq!(row.left_only, 0);
        assert_eq!(row.right_only, 1);
        assert_eq!(row.jaccard, 0.0);

        assert!(db.agreement_rows("missing-run")?.is_empty());

        fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn reports_ignore_running_runs() -> Result<()> {
        let path = temp_db_path("agreement-running-run");
        let db = BenchmarkDb::open(&path)?;
        insert_run(&db, "run-1", "2026-05-08T00:00:00Z", "running")?;
        insert_dataset_subset(&db, "run-1", "dataset-a")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "pairasm")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "fastp")?;
        insert_merged_read_record(&db, "run-1", "dataset-a", "pairasm", "read-1")?;

        assert!(db.agreement_rows("run-1")?.is_empty());
        assert!(db.tool_result_rows("run-1")?.is_empty());

        fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn records_merged_read_records_from_fastq() -> Result<()> {
        let path = temp_db_path("merged-read-records");
        let fastq_path = temp_fastq_path("merged-read-records");
        fs::write(
            &fastq_path,
            "@read/1 extra metadata\nACGT\n+\nIIII\n@other/2\nAC\n+\nII\n",
        )?;
        let db = BenchmarkDb::open(&path)?;
        insert_run(&db, "run-1", "2026-05-08T00:00:00Z", "completed")?;
        insert_dataset_subset(&db, "run-1", "dataset-a")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "pairasm")?;
        let subset = SubsetMetadata {
            name: "dataset-a".to_owned(),
            accession: "DRR000000".to_owned(),
            read_pairs: 2,
            r1: PathBuf::from("/tmp/r1.fastq"),
            r2: PathBuf::from("/tmp/r2.fastq"),
        };

        let inserted =
            db.insert_merged_read_records_from_fastq("run-1", &subset, Tool::Pairasm, &fastq_path)?;

        assert_eq!(inserted, 2);
        assert_eq!(
            merged_read_record_count(&db, "run-1", "dataset-a", "pairasm")?,
            2
        );
        assert_eq!(
            merged_read_record_field(&db, "run-1", "dataset-a", "pairasm", "read", "merged_len")?,
            "4"
        );
        assert_eq!(
            merged_read_record_field(
                &db,
                "run-1",
                "dataset-a",
                "pairasm",
                "read",
                "output_header"
            )?,
            "@read/1 extra metadata"
        );

        fs::remove_file(path)?;
        fs::remove_file(fastq_path)?;
        Ok(())
    }

    #[test]
    fn failed_tool_evidence_replacement_rolls_back_prior_rows() -> Result<()> {
        let path = temp_db_path("replace-tool-evidence-rollback");
        let fastq_path = temp_fastq_path("replace-tool-evidence-rollback");
        fs::write(&fastq_path, "@new-read/1\nA\n+\n!\n")?;
        let db = BenchmarkDb::open(&path)?;
        insert_run(&db, "run-1", "2026-05-08T00:00:00Z", "completed")?;
        insert_dataset_subset(&db, "run-1", "dataset-a")?;
        insert_tool_execution(&db, "run-1", "dataset-a", "pairasm")?;
        insert_artifact(&db, "run-1", "dataset-a", "pairasm", "command")?;
        insert_merged_read_record(&db, "run-1", "dataset-a", "pairasm", "old-read")?;
        let subset = SubsetMetadata {
            name: "dataset-a".to_owned(),
            accession: "DRR000000".to_owned(),
            read_pairs: 2,
            r1: PathBuf::from("/tmp/r1.fastq"),
            r2: PathBuf::from("/tmp/r2.fastq"),
        };
        let timing = TimingRecord {
            mean_s: 2.0,
            median_s: 2.0,
            stddev_s: None,
            min_s: 2.0,
            max_s: 2.0,
            user_s: 1.0,
            system_s: 1.0,
        };
        let execution = ToolExecutionRecord {
            run_key: "run-1",
            subset: &subset,
            tool: Tool::Pairasm,
            command_string: "new-command",
            execution_order: 1,
            timing: &timing,
            merged_reads: 2,
            expected_merged_read_records: 2,
            output_dir: Path::new("/tmp/new-out"),
        };
        let artifacts = [ArtifactRecord {
            kind: ArtifactKind::Command,
            path: Path::new("/tmp/new-command.sh"),
            requirement: ArtifactRequirement::Required,
        }];

        let error = match db.replace_tool_evidence(&execution, &artifacts, &fastq_path) {
            Ok(_) => panic!("count mismatch should fail"),
            Err(error) => error,
        };

        assert!(
            error.to_string().contains("merged read count changed"),
            "unexpected error: {error:?}"
        );
        assert_eq!(
            tool_execution_count(&db, "run-1", "dataset-a", "pairasm")?,
            1
        );
        assert_eq!(artifact_count(&db, "run-1", "dataset-a", "pairasm")?, 1);
        assert_eq!(
            merged_read_record_count(&db, "run-1", "dataset-a", "pairasm")?,
            1
        );
        assert_eq!(
            merged_read_record_field(&db, "run-1", "dataset-a", "pairasm", "old-read", "read_id")?,
            "old-read"
        );

        fs::remove_file(path)?;
        fs::remove_file(fastq_path)?;
        Ok(())
    }

    #[test]
    fn open_existing_rejects_missing_database_without_creating_it() {
        let path = temp_db_path("missing-report-db");

        let error = match BenchmarkDb::open_existing(&path) {
            Ok(_) => panic!("missing DB should fail"),
            Err(error) => error,
        };

        assert!(
            error
                .to_string()
                .contains("benchmark database does not exist")
        );
        assert!(!path.exists());
    }

    fn temp_db_path(test_name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "pairasm-benches-{test_name}-{}.duckdb",
            uuid::Uuid::new_v4()
        ))
    }

    fn temp_fastq_path(test_name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "pairasm-benches-{test_name}-{}.fastq",
            uuid::Uuid::new_v4()
        ))
    }

    fn insert_run(db: &BenchmarkDb, run_key: &str, created_at: &str, status: &str) -> Result<()> {
        let completed_at = (status == "completed").then_some(created_at);
        db.connection.execute(
            "INSERT INTO benchmark_runs (
                run_key, created_at, status, completed_at, benchmark_mode, read_pairs, replicates,
                threads, run_dir, config_path, data_root
            ) VALUES (?1, ?2, ?3, ?4, 'default-user', 1, 1, 1, '/tmp/run', '/tmp/config', '/tmp/data')",
            params![run_key, created_at, status, completed_at],
        )?;
        Ok(())
    }

    fn insert_dataset_subset(db: &BenchmarkDb, run_key: &str, dataset_name: &str) -> Result<()> {
        db.connection.execute(
            "INSERT INTO dataset_subsets (
                run_key, dataset_name, accession, read_pairs, r1_path, r2_path, r1_bytes, r2_bytes
            ) VALUES (?1, ?2, 'DRR000000', 1, '/tmp/r1.fastq', '/tmp/r2.fastq', 10, 10)",
            params![run_key, dataset_name],
        )?;
        Ok(())
    }

    fn insert_tool_execution(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
    ) -> Result<()> {
        db.connection.execute(
            "INSERT INTO tool_executions (
                run_key, dataset_name, tool, command, execution_order, hyperfine_mean_s, hyperfine_median_s,
                hyperfine_stddev_s, hyperfine_min_s, hyperfine_max_s, hyperfine_user_s,
                hyperfine_system_s, merged_reads, merged_read_records, output_dir
            ) VALUES (?1, ?2, ?3, ?3, 0, 1.0, 1.0, NULL, 1.0, 1.0, 0.5, 0.5, 0, 0, '/tmp/out')",
            params![run_key, dataset_name, tool],
        )?;
        Ok(())
    }

    fn insert_artifact(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
        artifact_kind: &str,
    ) -> Result<()> {
        db.connection.execute(
            "INSERT INTO artifacts (run_key, dataset_name, tool, artifact_kind, path, required, bytes)
             VALUES (?1, ?2, ?3, ?4, '/tmp/old-artifact', TRUE, NULL)",
            params![run_key, dataset_name, tool, artifact_kind],
        )?;
        Ok(())
    }

    fn insert_merged_read_record(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
        read_id: &str,
    ) -> Result<()> {
        db.connection.execute(
            "INSERT INTO merged_read_records (
                run_key, dataset_name, tool, read_id, output_header, merged_len,
                avg_qual, min_qual, max_qual, sequence_hash, quality_hash
            ) VALUES (?1, ?2, ?3, ?4, ?4, 10, 40.0, 40, 40, 'seq-hash', 'qual-hash')",
            params![run_key, dataset_name, tool, read_id],
        )?;
        Ok(())
    }

    fn merged_read_record_count(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
    ) -> Result<i64> {
        Ok(db.connection.query_row(
            "SELECT COUNT(*) FROM merged_read_records WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
            params![run_key, dataset_name, tool],
            |row| row.get(0),
        )?)
    }

    fn tool_execution_count(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
    ) -> Result<i64> {
        Ok(db.connection.query_row(
            "SELECT COUNT(*) FROM tool_executions WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
            params![run_key, dataset_name, tool],
            |row| row.get(0),
        )?)
    }

    fn artifact_count(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
    ) -> Result<i64> {
        Ok(db.connection.query_row(
            "SELECT COUNT(*) FROM artifacts WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3",
            params![run_key, dataset_name, tool],
            |row| row.get(0),
        )?)
    }

    fn merged_read_record_field(
        db: &BenchmarkDb,
        run_key: &str,
        dataset_name: &str,
        tool: &str,
        read_id: &str,
        field: &str,
    ) -> Result<String> {
        let sql = match field {
            "read_id" => {
                "SELECT read_id FROM merged_read_records WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3 AND read_id = ?4"
            },
            "merged_len" => {
                "SELECT merged_len::TEXT FROM merged_read_records WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3 AND read_id = ?4"
            },
            "output_header" => {
                "SELECT output_header FROM merged_read_records WHERE run_key = ?1 AND dataset_name = ?2 AND tool = ?3 AND read_id = ?4"
            },
            _ => panic!("unsupported merged read record field: {field}"),
        };
        Ok(db
            .connection
            .query_row(sql, params![run_key, dataset_name, tool, read_id], |row| {
                row.get(0)
            })?)
    }
}
