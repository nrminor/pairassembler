use thiserror::Error;

/// Custom Result type for libpairassembly operations, wrapping the custom [`enum@Error`] type
#[allow(clippy::absolute_paths)]
pub type Result<T> = std::result::Result<T, Error>;

/// The main error type for `libpairassembly`, encompassing all possible error cases
/// that can occur during binary sequence operations.
#[allow(clippy::enum_variant_names)]
#[derive(Error, Debug)]
#[error(transparent)]
pub enum Error {
    /// Errors related to converting external record types, e.g., the noodles FASTQ record type, into
    /// an internal representation that can be used for pairing, overlapping, etc.
    #[error(
        "Error converting external sequence read types into internal data structures for pairing:

{0}

`libpairassembly` currently uses a closed type system, meaning that users cannot themselves add \
support for additional FASTQ record types. If you would like support to be made open, please open \
an issue at https://github.com/nrminor/pairassembler/issues or submit a PR!"
    )]
    ConversionError(#[from] ConversionError),

    /// Errors related to input or output of data, which is to say errors at the earliest or latest stages of
    /// data processing with `libpairassembly`.
    #[error(
        "Error encountered during an attempt to input or output data from the `libpairassembly` workflow: {0}"
    )]
    InputOutputError(#[from] InputOutputError),

    /// Errors related to pairing two sequence reads as mates, which must occur successfully before
    /// overlapping, validation, merging, and correction can occur.
    #[error("Error encountered when attempting to pair sequence reads as mates: {0}")]
    PairingError(#[from] PairingError),

    /// Errors related to overlapping two sequence read mates.
    #[error("Error encountered when attempting to find overlapping bases between mates: {0}")]
    OverlapError(#[from] OverlapError),

    /// Errors related to validating potential overlap between sequence read mates.
    #[error("Error encountered when validating a successful overlap between two read mates: {0}")]
    ValidationError(#[from] ValidationError),

    /// Errors related to merging a pair of reads.
    #[error("Error encountered while merging a pair of reads.")]
    MergeError(#[from] MergeError),

    /// Errors related to overlap correction based on available information in the mates' quality scores.
    #[error(
        "Error encountered while using base-calls and quality scores in a validated overlap to \
        correct quality scores:

        {0}"
    )]
    CorrectionError(#[from] CorrectionError),

    /// Generic errors for other unexpected situations handled with anyhow
    #[error("Generic error: {0}")]
    AnyhowError(#[from] anyhow::Error),

    /// Generic errors for other unexpected situations handled with ColorEyre
    #[error("Generic error: {0}")]
    ColorEyreError(#[from] color_eyre::Report),
}

#[derive(Debug, Error)]
pub enum ConversionError {
    #[error("Failed to construct output record type from corrected parts: {0}")]
    RecordConstruction(String),
}
pub use ConversionError::*;

#[derive(Debug, Error)]
pub enum InputOutputError {
    #[error(
        "Mismatched sequence length and quality score length encountered:

{1} bases:     {0}
{3} qualities: {2}

FASTQ format requires that one Phred quality score is provided per base. Without this information, \
`libpairassembly` cannot proceed.
"
    )]
    SequenceQualityLengthMismatch(String, usize, String, usize),
}
pub use InputOutputError::*;

#[derive(Debug, Error)]
pub enum PairingError {
    /// Invalid FASTQ ID entry encountered.
    #[error(
        r"Invalid ID encountered that cannot be used to identify paired read mates:

{0}

Please make sure that sequence reads in the input FASTQ match the following template:

```
>{{id_prefix}}{{id}}{{id_suffix}} {{comment_prefix}}{{comment}}{{comment_suffix}}
{{sequence}}
+{{plus_prefix}}{{plus}}{{plus_suffix}}
{{quality}}
```
"
    )]
    InvalidId(String), // NOTE: Because this error type contains a string, it should be constructed lazily with something like `.unwrap_or_else(|error| error.to_string())`

    /// Error for when an attempt is made to incorrectly pair reads with different ID information.
    #[error(
        "Attempted to pair reads with unmatched IDs: {0} and {1} are not mates, and thus there is no \
good reason to expect that their bases will overlap or that their quality scores will meaningfully \
reflect the same template molecule."
    )]
    UnmatchedIds(String, String), // Same as above with allocations.

    /// Error for when a read is paired with itself
    #[error(
        "Attempt made to recursively pair a read with itself, a self-referential mate. Its header is '{0}'."
    )]
    RecursivePairing(String),
}
pub use PairingError::*;

#[derive(Debug, Error)]
pub enum OverlapError {
    /// Expected runtime outcome for pairs that do not overlap under current overlap settings.
    #[error("No overlap found for paired reads under current overlap settings.")]
    NoOverlapFound,

    /// Error for when an overlap is found that is below a minimum length. In most cases, this need not
    /// be escalated to the level of error. Instead, overlaps should be wrapped in an Option, where a
    /// failure to find an overlap is simply None, as a pair of reads without an overlap is a normal,
    /// expected outcome with paired-end sequencing.
    #[error("Overlap length {found} below minimum required {required}.")]
    OverlapBelowMinimum { found: usize, required: usize },

    /// Error for when an invalid base character, which is to say a character that is not supported
    /// in the FASTQ format specification, slipped through, leading to a reverse complement that is
    /// shorter than the original template.
    #[error("Reverse complement length mismatch: original = {original}, revcomp = {revcomp}")]
    ReverseComplementLengthMismatch { original: usize, revcomp: usize },

    /// Error for when the overlap algorithm has erroneously begun comparing bases outside the bounds
    /// of a read.
    #[error("Index {index} out of bounds for {read} (length {length}).")]
    IndexOutOfBounds {
        read: &'static str,
        index: usize,
        length: usize,
    },

    /// Error for when an overlap-oriented mate has different sequence and quality lengths.
    #[error(
        "sequence/quality length mismatch in overlap-oriented mate '{mate}': sequence={seq_len}, quality={qual_len}"
    )]
    OrientedPairSequenceQualityLengthMismatch {
        mate: &'static str,
        seq_len: usize,
        qual_len: usize,
    },

    /// Error for when an invalid overlap length that is longer than either read or shorter than the required minimum has somehow slipped through the cracks.
    #[error(
        "Invalid overlap length: computed length {computed} with bounds read1 = {read1_len}, read2 \
= {read2_len}, and min required = {min_required}"
    )]
    InvalidOverlapLength {
        computed: usize,
        read1_len: usize,
        read2_len: usize,
        min_required: usize,
    },

    /// Error for when an overlap starting from the 5' end of the pair and an overlap starting from
    /// the 3' end of the pair are both found *and* where both overlaps have the same rate of mismatches,
    /// which is to say the same number of mismatches divided by overlap length. We expect this to be
    /// a very rare occurrence.
    #[error(
        "Pair with ambiguous overlap status caused by two overlaps of equivalent quality \
(mismatches={diff}, overlap_len={overlap_len}). `libpairassembly` does not support such overlap \
ties at this time, though please raise an issue at https://github.com/nrminor/pairassembler/issues \
or submit a PR if this support is important to your use case."
    )]
    OverlapTie { diff: usize, overlap_len: usize },
}
pub use OverlapError::*;

#[derive(Debug, Error)]
pub enum ValidationError {
    /// Error for when an overlap was found, but is of insufficient length given the sequence
    /// complexity of the paired read mates.
    #[error(
        "The observed overlap length between two mated reads, {observed_overlap_len}, is insufficient \
 with the provided parameters. Overlaps with the provided K of {k} and minimum complexity score of {min_complexity_score} \
 must be at least {min_overlap_len} bases. As such, this overlap can justifiably be excluded from merging."
    )]
    InsufficientOverlapLength {
        observed_overlap_len: usize,
        min_overlap_len: usize,
        min_complexity_score: usize,
        k: usize,
    },

    /// Error for when a given overlap is assessed in terms of its sequence complexity as well as
    /// its rate of mismatches over the overlap's length and found to have too many mismatches.
    #[error(
        "Overlap encountered with a mismatch rate, {observed_error_rate}, that was higher than the \
 maximum expected error rate for the overlap, {maximum_expected_error_rate}, given its sequence \
 complexity, the provided K of {k}, and the provided minimum complexity score {min_complexity_score}. As such, this \
 overlap can justifiably be excluded from merging."
    )]
    ExcessiveObservedMismatchRate {
        min_complexity_score: usize,
        k: usize,
        observed_error_rate: f32,
        maximum_expected_error_rate: f32,
    },
}
pub use ValidationError::*;

#[derive(Debug, Error)]
pub enum MergeError {
    /// Error for when forward and reverse overlap windows have incompatible lengths.
    #[error("merge overlap windows length mismatch: fwd={fwd_len}, rev={rev_len}")]
    OverlapWindowLengthMismatch { fwd_len: usize, rev_len: usize },

    /// Error for when any merge section has mismatched sequence and quality lengths.
    #[error(
        "sequence/quality length mismatch in merge section '{section}': sequence={seq_len}, quality={qual_len}"
    )]
    MergeSequenceQualityLengthMismatch {
        section: &'static str,
        seq_len: usize,
        qual_len: usize,
    },

    /// Error for when merge receives an empty overlap window.
    #[error("overlap length must be greater than zero for merge")]
    EmptyOverlapWindow,

    /// Error for when the final merged read is length expected when summing the left overhang, the
    /// right overhang, and the overlap in between the two.
    #[error("Total merged read length ({actual}) does not match computed length ({expected}).")]
    MergedLengthMismatch { expected: usize, actual: usize },

    /// Error for when overlap provenance cannot satisfy expected overlap lengths.
    #[error(
        "merge provenance overlap length ({overlap_len}) does not match provenance vectors forward={fwd_len}, reverse={rev_len}"
    )]
    ProvenanceLengthMismatch {
        overlap_len: usize,
        fwd_len: usize,
        rev_len: usize,
    },

    /// Error for when merge policy rejects an equal-quality overlap base disagreement.
    #[error(
        "equal-quality overlap disagreement rejected at offset {offset}: forward={fwd_base:?}, reverse={rev_base:?}, quality={quality}"
    )]
    EqualQualityBaseDisagreement {
        offset: usize,
        fwd_base: u8,
        rev_base: u8,
        quality: u8,
    },
}
pub use MergeError::*;

#[derive(Debug, Error)]
pub enum CorrectionError {
    /// Error for when the consensus sequence and quality slices are of different lengths, which
    /// will make it impossible to run quality score correction for all bases.
    #[error("Mismatch between consensus sequence and quality vector: {seq_len} != {qual_len}.")]
    ConsensusLengthMismatch { seq_len: usize, qual_len: usize },

    /// Error for catching an invalid base.
    #[error("Invalid base encountered during correction.")]
    InvalidBase(u8),

    // TODO: Support for Phred+64 ASCII range?
    /// Error for catching an invalid quality score in the Phred+33 ASCII range.
    #[error("Quality score out of Phred+33 ASCII range: {0}")]
    InvalidQualityScore(u8),

    /// Error for when a quality score computation over/underflows the range of valid Phred scores.
    #[error("Overflow or underflow during Phred score computation.")]
    QualityScoreComputationError,

    /// Error for when slices of bases and slices of quality scores are not the same length.
    #[error("Unexpected mismatch in number of base-pairs during correction.")]
    AlignmentLengthMismatch { seq_len: usize, qual_len: usize },
}
pub use CorrectionError::*;
