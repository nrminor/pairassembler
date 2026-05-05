use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::model::Tool;

const DEFAULT_CONFIG: &str = "benches/config/datasets.tsv";
const DEFAULT_DATA_ROOT: &str = "benches/data";
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
    /// Check external benchmark tools and print versions.
    Check,
    /// Fetch configured paired FASTQs from ENA.
    Fetch(CommonOptions),
    /// Prepare deterministic first-N-pair FASTQ subsets.
    Prepare(PrepareOptions),
    /// Run pairasm and competitor merge tools through hyperfine.
    Run(RunOptions),
    /// Summarize hyperfine run artifacts into a TSV table.
    Summarize(SummarizeOptions),
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
    #[arg(long, default_value_t = OutputCompression::Plain)]
    pub output_compression: OutputCompression,
    #[arg(long, default_value_t = BenchmarkMode::DefaultUser)]
    pub mode: BenchmarkMode,
}

#[derive(Debug, Clone, Parser)]
pub struct SummarizeOptions {
    #[arg(long, default_value = DEFAULT_RUNS_ROOT)]
    pub runs_root: PathBuf,
    #[arg(long)]
    pub run_dir: Option<PathBuf>,
    #[arg(long, default_value_t = true)]
    pub latest: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum BenchmarkMode {
    /// Minimal-thought CLI defaults from paired R1/R2 FASTQs.
    DefaultUser,
    /// Explicitly tuned/demo settings that make tool policies more comparable.
    TunedComparability,
}

impl BenchmarkMode {
    pub fn name(self) -> &'static str {
        match self {
            Self::DefaultUser => "default-user",
            Self::TunedComparability => "tuned-comparability",
        }
    }
}

impl std::fmt::Display for BenchmarkMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputCompression {
    Plain,
    Gzip,
}

impl std::fmt::Display for OutputCompression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Plain => f.write_str("plain"),
            Self::Gzip => f.write_str("gzip"),
        }
    }
}
