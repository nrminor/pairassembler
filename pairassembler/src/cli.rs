use clap::{
    Parser,
    builder::{
        Styles,
        styling::{AnsiColor, Effects},
    },
};

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
    pub verbose: clap_verbosity_flag::Verbosity,

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
    #[arg(long, default_value_t = 2, help_heading = "Overlap Settings")]
    pub overlap_diff_max: usize,

    /// Minimum number of bases that must overlap before an overlap can be accepted.
    #[arg(long, default_value_t = 30, help_heading = "Overlap Settings")]
    pub min_overlap: usize,

    /// Maximum mismatch fraction allowed for longer overlaps.
    #[arg(long, default_value_t = 0.2, help_heading = "Overlap Settings")]
    pub diff_percent_max: f32,

    /// Minimum base comparisons required before overlap thresholds are meaningful.
    #[arg(long, default_value_t = 50, help_heading = "Overlap Settings")]
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
        default_value_t = 39,
        help_heading = "Validation Settings"
    )]
    pub min_complexity_score: usize,

    /// Merge reads without overlap-based quality correction.
    #[arg(
        short = 'n',
        long,
        default_value_t = false,
        help_heading = "Correction Settings"
    )]
    pub no_correct: bool,
}

// Configure Clap help menu colors.
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Yellow.on_default());
