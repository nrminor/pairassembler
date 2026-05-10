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
//!
//! Most users should start with [`Assembler`] or the common exports in [`prelude`].
//!
//! # Quick start
//!
//! ```rust
//! use libpairassembly::prelude::*;
//!
//! # fn main() -> libpairassembly::Result<()> {
//! let forward = SequenceRead::try_new(
//!     "read-1",
//!     "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
//!     "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
//! )?;
//! let reverse = SequenceRead::try_new(
//!     "read-1",
//!     "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
//!     "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
//! )?;
//!
//! let pair = PairInput::new(forward, reverse);
//! let mut assembler = Assembler::default();
//!
//! let merged = assembler
//!     .process_pair(&pair)?
//!     .expect("this fixture is a full-length overlap");
//!
//! assert_eq!(merged.id(), "read-1");
//! assert_eq!(merged.sequence_bytes().len(), merged.quality_bytes().len());
//! # Ok(())
//! # }
//! ```
//!
//! Pairs without an acceptable overlap return `Ok(None)`, not `Err`:
//!
//! ```rust
//! use libpairassembly::prelude::*;
//!
//! # fn main() -> libpairassembly::Result<()> {
//! let pair = PairInput::new(
//!     SequenceRead::try_new(
//!         "read-1",
//!         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
//!         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
//!     )?,
//!     SequenceRead::try_new(
//!         "read-1",
//!         "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
//!         "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
//!     )?,
//! );
//!
//! let mut assembler = Assembler::default();
//! assert!(assembler.process_pair(&pair)?.is_none());
//! # Ok(())
//! # }
//! ```
mod overlap;

mod validate;

mod merge;

mod correct;

pub mod c_abi;

pub mod read;

pub mod assembler;

pub mod errors;
pub use errors::{Error, Result};

pub mod prelude;

#[cfg(test)]
pub(crate) mod test_fixtures;

pub use assembler::*;
pub use prelude::*;
