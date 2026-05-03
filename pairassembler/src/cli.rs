use clap::{
    Parser, Subcommand,
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
like those produced by some Illumina platforms. It can use this information to merge read mates into \
consensus reads and also to correct quality scores where both reads overlap, simplifying and improving \
downstream genomic analysis in the process. Reads can be input in FASTQ, gzip-compressed FASTQ, and \
BINSEQ format and can be output in the same formats in addition to FASTA and gzip-compressed FASTA. \
";

#[derive(Parser)]
#[clap(name = "pairasm", about = INFO, version)]
#[command(styles = STYLES)]
pub struct Cli {
    #[command(flatten)]
    pub verbose: clap_verbosity_flag::Verbosity,

    #[command(subcommand)]
    pub command: Option<Commands>,

    /// Crash early if read mates aren't in the same order
    #[arg(long)]
    pub strict: bool,
}

#[derive(Subcommand)]
pub enum Commands {
    #[clap(
        about = "Identify overlaps between read mates in input files and merge them into a single read per pair. Unmerged reads can optionally be written out to a separate output file.",
        visible_aliases = &["m", "pair", "p", "pr"],
        aliases = &[ "mrge", "mrege", "piar"],
    )]
    Merge {
        /// First input file to draw potentially overlapping reads from. If only a single input file
        /// is provided, mates are assumed to be interleaved.
        #[arg(short = '1', long, required = true, help_heading = "Inputs")]
        input1: String,

        /// Optional second input file to draw potentially overlapping reads from. If not provided,
        /// reads in the first input file will assumed to be interleaved mates.
        #[arg(short = '2', long, required = false, help_heading = "Inputs")]
        input2: Option<String>,

        /// Output file to write merged reads into. If not specified, defaults to standard output.
        /// File extension is used to determine format and compression codec.
        #[arg(short, long, required = false, help_heading = "Outputs")]
        output_file: Option<String>,

        /// Optional output file to write unmerged reads into.
        #[arg(short, long, required = false, help_heading = "Outputs")]
        unmerged_out: Option<String>,

        /// The maximum number of basecall mismatches permitted before a potential overlap between mates is rejected.
        #[arg(
            long,
            required = false,
            default_value_t = 2,
            help_heading = "Overlap Settings"
        )]
        overlap_diff_max: usize,

        /// The minimum number of bases that must overlap for an overlap to be accepted.
        #[arg(
            long,
            required = false,
            default_value_t = 30,
            help_heading = "Overlap Settings"
        )]
        min_overlap: usize,

        /// Complementary with `--overlap-diff-max`; the percentage of a very long overlap that is
        /// allowed to be mismatches before that overlap is rejected.
        #[arg(
            long,
            required = false,
            default_value_t = 0.2,
            help_heading = "Overlap Settings"
        )]
        diff_percent_max: f32,

        /// The minimum number of base calls that must be compared before an overlap between mates can
        /// be accepted.
        #[arg(
            long,
            required = false,
            default_value_t = 50,
            help_heading = "Overlap Settings"
        )]
        min_comparisons: usize,

        /// Kmer length to use when assessing the information content of an overlap, defaulting to
        /// 3. This value should generally not be changed; increasing it will allow more
        /// overlaps to pass validation.
        #[arg(
            short,
            long = "kmer-size",
            required = false,
            default_value_t = 3,
            help_heading = "Validation Settings"
        )]
        k: usize,

        /// Minimum k-mer complexity score used when determining whether a detected overlap is
        /// informative enough to trust. 39 will be appropriate in most cases; lower values allow
        /// more overlaps and higher values allow fewer.
        #[arg(
            short = 'E',
            long,
            required = false,
            default_value_t = 39,
            help_heading = "Validation Settings"
        )]
        min_complexity_score: usize,

        /// Turn off quality score correction
        #[arg(
            short,
            long,
            required = false,
            default_value_t = false,
            help_heading = "Correction Settings"
        )]
        no_correct: bool,
    },

    #[clap(
        about = "Identify overlaps between read mates and use overlapping base-calls to correct quality scores on each mate without merging them. Corrected reads will be interleaved into a single output file.",
        visible_aliases = &["c", "fix", "adjust", "a"],
        aliases = &["corct", "corrcet", "fxi"],
    )]
    Correct {
        /// First input file to draw potentially overlapping reads from. If only a single input file
        /// is provided, mates are assumed to be interleaved.
        #[arg(short = '1', long, required = true, help_heading = "Inputs")]
        input1: String,

        /// Optional second input file to draw potentially overlapping reads from. If not provided,
        /// reads in the first input file will assumed to be interleaved mates.
        #[arg(short = '2', long, required = false, help_heading = "Inputs")]
        input2: Option<String>,

        /// Output file to write merged reads into. If not specified, defaults to standard output.
        /// File extension is used to determine format and compression codec.
        #[arg(short, long, required = false, help_heading = "Outputs")]
        output_file: Option<String>,

        /// Optional output file to write unmerged reads into.
        #[arg(short, long, required = false, help_heading = "Outputs")]
        unmerged_out: Option<String>,

        /// The maximum number of basecall mismatches permitted before a potential overlap between
        /// mates is rejected.
        #[arg(
            long,
            required = false,
            default_value_t = 2,
            help_heading = "Overlap Settings"
        )]
        overlap_diff_max: usize,

        /// The minimum number of bases that must overlap for an overlap to be accepted.
        #[arg(
            long,
            required = false,
            default_value_t = 30,
            help_heading = "Overlap Settings"
        )]
        min_overlap: usize,

        /// Complementary with `--overlap-diff-max`; the percentage of a very long overlap that is
        /// allowed to be mismatches before that overlap is rejected.
        #[arg(
            long,
            required = false,
            default_value_t = 0.2,
            help_heading = "Overlap Settings"
        )]
        diff_percent_max: f32,

        /// The minimum number of base calls that must be compared before an overlap between mates can
        /// be accepted.
        #[arg(
            long,
            required = false,
            default_value_t = 50,
            help_heading = "Overlap Settings"
        )]
        min_comparisons: usize,

        /// Kmer length to use when assessing the information content of an overlap, defaulting to
        /// 3. This value should generally not be changed; increasing it will allow more
        /// overlaps to pass validation.
        #[arg(
            short,
            long = "kmer_size",
            required = false,
            default_value_t = 3,
            help_heading = "Validation Settings"
        )]
        k: usize,

        /// Minimum k-mer complexity score used when determining whether a detected overlap is
        /// informative enough to trust. 39 will be appropriate in most cases; lower values allow
        /// more overlaps and higher values allow fewer.
        #[arg(
            short = 'E',
            long,
            required = false,
            default_value_t = 39,
            help_heading = "Validation Settings"
        )]
        min_complexity_score: usize,
    },

    #[clap(
        about = "NOT YET IMPLEMENTED. Check that sequencing read mates are in the correct order in the input files. This subcommand is intended to be run prior to other `pairasm` subcommands to prevent wasted computation. This subcommand will also identify orphans that may disrupt pairing.",
        visible_aliases = &["v", "val", "check", "ch"],
    )]
    Validate {
        /// First input file to draw potentially overlapping reads from. If only a single input file
        /// is provided, mates are assumed to be interleaved.
        #[arg(short = '1', long, required = true, help_heading = "Inputs")]
        input1: String,

        /// Optional second input file to draw potentially overlapping reads from. If not provided,
        /// reads in the first input file will assumed to be interleaved mates.
        #[arg(short = '2', long, required = false, help_heading = "Inputs")]
        input2: Option<String>,
    },

    #[clap(
        about = "NOT YET IMPLEMENTED. Validate that sequence read mates are in the correct order in the input files and sort them if not.",
        visible_aliases = &["so", "s", "ord", "order", "reorder", "repair"],
    )]
    Sort {
        /// First input file to draw potentially overlapping reads from. If only a single input file
        /// is provided, mates are assumed to be interleaved.
        #[arg(short = '1', long, required = true, help_heading = "Inputs")]
        input1: String,

        /// Optional second input file to draw potentially overlapping reads from. If not provided,
        /// reads in the first input file will assumed to be interleaved mates.
        #[arg(short = '2', long, required = false, help_heading = "Inputs")]
        input2: Option<String>,
    },
}

// Configures Clap v3-style help menu colors
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Yellow.on_default());
