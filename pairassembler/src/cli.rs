use std::{io, io::IsTerminal, path::PathBuf};

use clap::{
    ArgAction, Parser,
    builder::{
        Styles,
        styling::{AnsiColor, Effects},
    },
};
use color_eyre::eyre::{Result, eyre};
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use crate::progress::ProgressMode;

pub const INFO: &str = "

▄▄▄▄  ▗▞▀▜▌▄  ▄▄▄ ▗▞▀▜▌ ▄▄▄  ▄▄▄ ▗▞▀▚▖▄▄▄▄  ▗▖   █ ▗▞▀▚▖ ▄▄▄ 
█   █ ▝▚▄▟▌▄ █    ▝▚▄▟▌▀▄▄  ▀▄▄  ▐▛▀▀▘█ █ █ ▐▌   █ ▐▛▀▀▘█    
█▄▄▄▀      █ █         ▄▄▄▀ ▄▄▄▀ ▝▚▄▄▖█   █ ▐▛▀▚▖█ ▝▚▄▄▖█    
█                                           ▐▙▄▞▘█           
▀
pairassembler (v0.1.0)
------------------------------------------------------------
PairAssembler, called with `pairasm` in the command line, identifies overlaps between paired reads \
like those produced by some Illumina platforms. It uses that information to merge read mates into \
consensus reads and optionally correct quality scores where both reads overlap.\
";

const EXAMPLES: &str = "Examples:
  Merge paired FASTQs and write merged reads:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o merged.fastq.gz

  Keep unmerged pairs in a separate file:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o merged.fastq.gz --unmerged-out unmerged.fastq.gz

  Tune overlap validation for more permissive merging:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --min-overlap 20 --min-complexity-score 30

  Merge detected overlaps without validation:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --no-validate

  Write a JSON run summary:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --summary run-summary.json

  Merge without overlap-based quality correction:
    pairasm -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --no-correct
";

#[derive(Parser)]
#[command(
    name = "pairasm",
    version,
    about = "Merge overlapping paired FASTQ reads.",
    long_about = INFO,
    after_help = EXAMPLES,
    styles = STYLES,
    override_usage = "pairasm -1 <R1.fastq[.gz]> -2 <R2.fastq[.gz]> [OPTIONS]"
)]
pub struct Cli {
    #[command(flatten)]
    pub verbosity: Verbosity,

    /// First FASTQ file containing the forward/R1 mates.
    #[arg(
        short = '1',
        long,
        visible_alias = "in1",
        value_name = "FASTQ",
        help_heading = "Inputs"
    )]
    pub input1: String,

    /// Second FASTQ file containing the reverse/R2 mates.
    #[arg(
        short = '2',
        long,
        visible_alias = "in2",
        value_name = "FASTQ",
        help_heading = "Inputs"
    )]
    pub input2: String,

    /// Output file for merged reads. If omitted, merged reads are written to standard output.
    /// File extension is used to determine format and compression codec.
    #[arg(short = 'o', long, value_name = "FASTQ", help_heading = "Outputs")]
    pub output_file: Option<String>,

    /// Optional output file for pairs that did not produce an accepted overlap.
    #[arg(short = 'u', long, value_name = "FASTQ", help_heading = "Outputs")]
    pub unmerged_out: Option<String>,

    /// Maximum basecall mismatches permitted before a potential overlap is rejected.
    #[arg(long, default_value_t = 5, help_heading = "Overlap Settings")]
    pub overlap_diff_max: usize,

    /// Minimum number of bases that must overlap before an overlap can be accepted.
    #[arg(long, default_value_t = 30, help_heading = "Overlap Settings")]
    pub min_overlap: usize,

    /// Maximum mismatch fraction allowed for longer overlaps.
    #[arg(long, default_value_t = 0.2, help_heading = "Overlap Settings")]
    pub diff_percent_max: f32,

    /// Minimum base comparisons required before overlap thresholds are meaningful.
    #[arg(long, default_value_t = 30, help_heading = "Overlap Settings")]
    pub min_comparisons: usize,

    /// K-mer length used when assessing overlap informativeness.
    #[arg(
        short = 'k',
        long = "kmer-size",
        default_value_t = 3,
        help_heading = "Validation Settings"
    )]
    pub k: usize,

    /// Minimum k-mer complexity score required to trust a detected overlap.
    #[arg(
        short = 'E',
        long,
        default_value_t = 30,
        help_heading = "Validation Settings"
    )]
    pub min_complexity_score: usize,

    /// Merge detected overlaps without overlap informativeness validation.
    #[arg(long, default_value_t = false, help_heading = "Validation Settings")]
    pub no_validate: bool,

    /// Merge reads without overlap-based quality correction.
    #[arg(
        short = 'n',
        long,
        default_value_t = false,
        help_heading = "Correction Settings"
    )]
    pub no_correct: bool,

    /// Number of complete pairs between progress updates.
    #[arg(long, default_value_t = 100_000, help_heading = "Reporting")]
    pub progress_every: u64,

    /// Write a JSON run summary to this path.
    #[arg(long, value_name = "JSON", help_heading = "Reporting")]
    pub summary: Option<PathBuf>,

    /// Fail after this many mate ID/order mismatches.
    #[arg(long, default_value_t = 3, help_heading = "Input Contract")]
    pub max_mate_id_mismatches: u64,
}

/// Logging and reporting verbosity flags.
#[derive(Clone, Copy, Debug, Default, Parser)]
#[command(about = None, long_about = None)]
pub struct Verbosity {
    /// Increase tracing verbosity (`-v` = INFO, `-vv` = DEBUG, `-vvv` = TRACE).
    #[arg(short = 'v', long = "verbose", action = ArgAction::Count, conflicts_with = "quiet")]
    pub verbose: u8,

    /// Reduce logs and reporting (`-q` = ERROR only, `-qq` = no logs or summary, `-qqq` = fully silent).
    #[arg(short = 'q', long = "quiet", action = ArgAction::Count, conflicts_with = "verbose")]
    pub quiet: u8,
}

/// Rendering and logging policy derived from verbosity flags and terminal state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct UiPolicy {
    pub log_level: Option<LevelFilter>,
    pub show_summary: bool,
    pub progress_mode: ProgressMode,
}

impl Cli {
    /// Initialize stderr-backed tracing using CLI verbosity and any `RUST_LOG` override.
    ///
    /// # Errors
    ///
    /// Returns an error when the global tracing subscriber cannot be initialized.
    pub fn init_tracing(&self) -> Result<()> {
        let Some(level) = self.ui_policy().log_level else {
            return Ok(());
        };

        let filter = EnvFilter::builder()
            .with_default_directive(level.into())
            .from_env_lossy();

        tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(io::stderr)
            .try_init()
            .map_err(|error| eyre!("failed to initialize tracing subscriber: {error}"))
    }

    /// Derive logging, progress, and summary behavior from verbosity flags.
    #[must_use]
    pub fn ui_policy(&self) -> UiPolicy {
        let stderr_is_tty = io::stderr().is_terminal();

        let (log_level, show_summary, show_progress) = match self.verbosity.quiet {
            0 => (
                Some(match self.verbosity.verbose {
                    0 => LevelFilter::WARN,
                    1 => LevelFilter::INFO,
                    2 => LevelFilter::DEBUG,
                    _ => LevelFilter::TRACE,
                }),
                true,
                true,
            ),
            1 => (Some(LevelFilter::ERROR), true, true),
            2 => (None, false, true),
            _ => (None, false, false),
        };

        let progress_mode = if !show_progress {
            ProgressMode::Off
        } else if stderr_is_tty {
            ProgressMode::Live
        } else {
            ProgressMode::Plain
        };

        UiPolicy {
            log_level,
            show_summary,
            progress_mode,
        }
    }
}

// Configure Clap help menu colors.
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Yellow.on_default());

#[cfg(test)]
mod tests {
    use tracing::level_filters::LevelFilter;

    use super::{Cli, UiPolicy, Verbosity};
    use crate::progress::ProgressMode;

    fn policy(verbose: u8, quiet: u8) -> UiPolicy {
        let cli = Cli {
            verbosity: Verbosity { verbose, quiet },
            input1: "r1.fastq".to_owned(),
            input2: "r2.fastq".to_owned(),
            output_file: None,
            unmerged_out: None,
            overlap_diff_max: 5,
            min_overlap: 30,
            diff_percent_max: 0.2,
            min_comparisons: 30,
            k: 3,
            min_complexity_score: 30,
            no_validate: false,
            no_correct: false,
            progress_every: 100_000,
            summary: None,
            max_mate_id_mismatches: 3,
        };
        cli.ui_policy()
    }

    #[test]
    fn ui_policy_defaults_to_warn_logs_and_reporting() {
        let policy = policy(0, 0);
        assert_eq!(policy.log_level, Some(LevelFilter::WARN));
        assert!(policy.show_summary);
    }

    #[test]
    fn ui_policy_maps_verbose_to_debug() {
        let policy = policy(2, 0);
        assert_eq!(policy.log_level, Some(LevelFilter::DEBUG));
        assert!(policy.show_summary);
    }

    #[test]
    fn ui_policy_maps_triple_quiet_to_silence() {
        let policy = policy(0, 3);
        assert_eq!(policy.log_level, None);
        assert!(!policy.show_summary);
        assert_eq!(policy.progress_mode, ProgressMode::Off);
    }
}
