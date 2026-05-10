use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::model::Tool;

const DEFAULT_CONFIG: &str = "benches/config/datasets.tsv";
const DEFAULT_DATA_ROOT: &str = "benches/data";
const DEFAULT_DB_PATH: &str = "benches/benchmarks.duckdb";
const DEFAULT_RUNS_ROOT: &str = "benches/runs";
const DEFAULT_READ_PAIRS: usize = 100_000;
const DEFAULT_REPLICATES: usize = 3;
const DEFAULT_THREADS: usize = 8;

#[derive(Debug, Parser)]
#[command(about = "Real-data comparative benchmarks for pairasm")]
pub struct Cli {
    #[command(subcommand)]
    pub command: BenchCommand,
}

#[derive(Debug, Subcommand)]
pub enum BenchCommand {
    /// Print a styled benchmark workflow phase banner.
    #[command(name = "workflow-phase", hide = true)]
    WorkflowPhase(WorkflowPhaseOptions),
    /// Check external benchmark tools and print versions.
    Check,
    /// Fetch configured paired FASTQs from ENA.
    Fetch(CommonOptions),
    /// Prepare deterministic first-N-pair FASTQ subsets.
    Prepare(PrepareOptions),
    /// Run pairasm and competitor merge tools through hyperfine.
    Run(RunOptions),
    /// Print benchmark reports from recorded results.
    Report(ReportOptions),
}

#[derive(Debug, Clone, Parser)]
pub struct WorkflowPhaseOptions {
    pub step: String,
    pub title: String,
}

#[derive(Debug, Clone, Parser)]
pub struct CommonOptions {
    #[arg(long, default_value = DEFAULT_CONFIG)]
    pub config: PathBuf,
    #[arg(long, default_value = DEFAULT_DATA_ROOT)]
    pub data_root: PathBuf,
}

#[derive(Debug, Clone, Parser)]
pub struct PrepareOptions {
    #[command(flatten)]
    pub common: CommonOptions,
    #[arg(long, default_value_t = DEFAULT_READ_PAIRS)]
    pub read_pairs: usize,
}

#[derive(Debug, Clone, Parser)]
pub struct RunOptions {
    #[command(flatten)]
    pub common: CommonOptions,
    #[arg(long, default_value = DEFAULT_RUNS_ROOT)]
    pub runs_root: PathBuf,
    #[arg(long, default_value = DEFAULT_DB_PATH)]
    pub db: PathBuf,
    #[arg(long, default_value_t = DEFAULT_READ_PAIRS)]
    pub read_pairs: usize,
    #[arg(long, default_value_t = DEFAULT_REPLICATES)]
    pub replicates: usize,
    #[arg(long, default_value_t = DEFAULT_THREADS)]
    pub threads: usize,
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "pairasm,fastp,bbmerge,vsearch"
    )]
    pub tools: Vec<Tool>,
    #[arg(long, default_value_t = BenchmarkMode::DefaultUser)]
    pub mode: BenchmarkMode,
}

#[derive(Debug, Clone, Parser)]
pub struct ReportOptions {
    #[command(subcommand)]
    pub command: ReportCommand,
}

#[derive(Debug, Clone, Subcommand)]
pub enum ReportCommand {
    /// Compare merged read-ID sets between tools.
    #[command(name = "read-id-overlap")]
    ReadIdOverlap(RunScopedReportOptions),
    /// Export per-tool timing/count results as TSV from DuckDB.
    #[command(name = "tool-results-tsv")]
    ToolResultsTsv(RunScopedReportOptions),
    /// Print per-tool timing/count results as a Markdown table from DuckDB.
    #[command(name = "timing-markdown")]
    TimingMarkdown(RunScopedReportOptions),
}

#[derive(Debug, Clone, Parser)]
pub struct RunScopedReportOptions {
    #[arg(long, default_value = DEFAULT_DB_PATH)]
    pub db: PathBuf,
    /// Run key to report. Defaults to the latest completed run for --mode.
    #[arg(long)]
    pub run: Option<String>,
    /// Benchmark mode to use when selecting the latest run.
    #[arg(long, default_value_t = BenchmarkMode::DefaultUser)]
    pub mode: BenchmarkMode,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum BenchmarkMode {
    /// Minimal-thought CLI defaults from paired R1/R2 FASTQs.
    DefaultUser,
    /// Explicitly tuned/demo settings that make tool policies more comparable.
    TunedComparability,
}

impl std::fmt::Display for BenchmarkMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DefaultUser => f.write_str("default-user"),
            Self::TunedComparability => f.write_str("tuned-comparability"),
        }
    }
}
