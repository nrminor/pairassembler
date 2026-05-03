//! `libpairassembly` is a Rust library crate for assembling and merging overlapping sequencing read
//! pairs into potentially higher quality consensus reads.
//!
//! The crate provides an [`assembler::Assembler`] API for pairing input records, finding acceptable
//! no-gap overlaps, optionally validating overlap informativeness, merging overlapping mates, and
//! applying overlap-aware quality-score correction. Pairs without an acceptable overlap are treated
//! as a normal biological outcome rather than an operational error.
//!
//! The main processing stages are organized across focused modules:
//!
//! 1. `read`: borrowed and owned read/pair domain types.
//! 2. `overlap`: overlap search parameters, oriented slices, and discovered overlap bounds.
//! 3. `validate`: overlap validation policies and retained validation metrics.
//! 4. `merge`: deterministic consensus construction from an overlap.
//! 5. `correct`: overlap-aware base and quality-score correction.
//! 6. `io`: optional FASTQ integration helpers behind feature flags.
//!
//! Most users should start with [`Assembler`] or the common exports in [`prelude`].
pub mod io;

mod overlap;

mod validate;

mod merge;

mod correct;

pub mod read;

pub mod assembler;

pub mod errors;
pub use errors::{Error, Result};

pub mod prelude;

#[cfg(test)]
pub(crate) mod test_fixtures;

pub use assembler::*;
pub use prelude::*;
