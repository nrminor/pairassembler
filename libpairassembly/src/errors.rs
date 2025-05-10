use thiserror::Error;

/// Custom Result type for libpairassembly operations, wrapping the custom [`Error`] type
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
        r#"Error converting external sequence read types into internal data structures for pairing:

{0}


`libpairassembly` currently uses a closed type system, meaning that users cannot themselves add support for additional FASTQ record types. If you would like support to be made open, please open an issue at https://github.com/nrminor/pairassembler/issues or submit a PR!"#
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
pub enum ConversionError {}
pub use ConversionError::*;

#[derive(Debug, Error)]
pub enum InputOutputError {
    #[error(r#"Mismatched sequence length and quality score length encountered:

{1} bases:     {0}
{3} qualities: {2}

FASTQ format requires that one Phred quality score is provided per base. Without this information, `libpairassembly` cannot proceed.
"#
    )]
    SequenceQualityLengthMismatch(String, usize, String, usize),
}
pub use InputOutputError::*;

#[derive(Debug, Error)]
pub enum PairingError {
    /// Invalid FASTQ ID entry encountered.
    #[error(
        r#"Invalid ID encountered that cannot be used to identify paired read mates:

{0}


Please make sure that sequence reads in the input FASTQ match the following template:

```
>{{id_prefix}}{{id}}{{id_suffix}} {{comment_prefix}}{{comment}}{{comment_suffix}}
{{sequence}}
+{{plus_prefix}}{{plus}}{{plus_suffix}}
{{quality}}
```
"#
    )]
    InvalidId(String), // NOTE: Because this error type contains a string, it should be constructed lazily with something like `.unwrap_or_else(|error| error.to_string())`

    /// Error for when an attempt is made to incorrectly pair reads with different ID information.
    #[error(
        "Attempted to pair reads with unmatched IDs: {0} and {1} are not mates, and thus there is not good reason to expect that their bases will overlap or that their quality scores will meaningfully reflect the same template molecule."
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

    /// Error for when an invalid overlap length that is longer than either read or shorter than the required minimum has somehow slipped through the cracks.
    #[error(
        "Invalid overlap length: computed length {computed} with bounds read1 = {read1_len}, read2 = {read2_len}, and min required = {min_required}"
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
        "Pair with ambiguous overlap status caused by two overlaps of equivalent quality, {0}, encountered. `libpairassembly` does not support such overlap ties at this type, though please raise an issue at https://github.com/nrminor/pairassembler/issues or submit a PR if this support is important to your use case."
    )]
    OverlapTie(f32),
}
pub use OverlapError::*;

#[derive(Debug, Error)]
pub enum ValidationError {}
pub use ValidationError::*;

#[derive(Debug, Error)]
pub enum MergeError {}
pub use MergeError::*;

#[derive(Debug, Error)]
pub enum CorrectionError {}
pub use CorrectionError::*;
