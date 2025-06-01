// dev allowances
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]
//
// crate-level lints
#![warn(
    // clippy::pedantic,
    clippy::perf,
    clippy::todo,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::complexity,
    clippy::correctness,
    clippy::absolute_paths,
    clippy::style
)]

//! `libpairassembly` is a Rust library crate for assembling and merging overlapping sequencing read
//! pairs (sometimes call mates) into potentially higher quality consensus reads. Particularly for
//! paired-end reads from Illumina instruments, pair assembly can both simplify and improve downstream
//! processing by reducing two files of read data ("R1" and "R2" FASTQ files) to a single file and by
//! improving quality scores when paired read basecalls agree where they overlap.
//!
//! `libpairassembly` takes inspiration from previous tools, including `fastp` for overlap discovery,
//! BBMerge for pre-merge overlap validation, and USEARCH/VSEARCH for quality score adjustment of
//! overlapping basecalls. This library aims to combine the best ideas from these three libraries
//! for the first time, exposing them as an idiomatic Rust API, a performant command line tool, and
//! an accessible, BioPython-compatible Python API. It also aims to streamline these tools, using
//! sensible defaults while only exposing parameters to tune behavior when they are likely to make
//! a difference in this particular use case (as opposed to providing a very large surface area
//! for tuning across large toolkits like VSEARCH or BBMap).
//!
//! In each of its interfaces, `libpairassembly` uses a standardized data flow for read-processing,
//! which is organized across modules:
//!
//! 1. `io`: Confirm that reads can be paired and find each read's mate.
//! 2. `validate`: Compute the information entropy for each pair to calibrate how large an overlap should be expected, a cool idea adopted from BBMerge.
//! 3. `overlap`: Find the longest overlap with the fewest mismatches.
//! 4. `correct`: Perform per-base Bayesian quality score correction as described in Edgar & Flyvbjerg (2015) and implemented in USEARCH and VSEARCH.
//! 5. `io` (again): Write merged reads into a new FASTQ, optionally writing unmerged reads to the same FASTQ or their own FASTQ.
//!
//! In addition to some of the methods from BBMerge, written in Java, and USEARCH/VSEARCH, written in C++,
//! `libpairassembly` brings with it a healthy dose of Rust goodness, including asynchronous read streaming, I/O that
//! is generic over a few FASTQ I/O crates, a few convenience macros, high performance, and a pipelined fluent interface
//! for library users.
//!

// TODO:
// - Even more docs, including doctests!
// - Benchmarks!
//

// generic record conversion interfaces with support for FASTQ I/O crates behind feature flags. The
// `libpairassembly` workflow starts and ends with symbols in the `io` module.
pub mod io;

// Like with data processing from particle colliders at CERN, `libpairassembly` works through a
// hierarchy of triage steps, each thinning the input data into cleaner, less noisy output data.
// The first step that performs the most data thinning involves the symbols from `overlap`, which
// provides functionality for finding the best overlaps between mated pairs. Reads that do not
// contain overlaps undergo no further processing, and are either done away with or written out.
mod overlap;

// Some overlaps are untrustworthy. The module `validate` uses as much available information from the
// reads as possible to determine whether found overlaps should be retained. Pairs with multiple
// equally valid overlaps, overlaps that are too short, or overlaps with too many mismatches are
// filtered out and henceforth treated as un-paired, joining reads that had no overlaps.
mod validate;

// Validate read overlaps are then merged in the `merge` module, with help from the `correct` module
// for adjusting quality scores (more on that below). A relatively simple method of finding overlaps
// is used, though local alignment with better support for indels may be added in the future.
mod merge;

// Sequencing reads with quality scores carry helpful information for read pair merging, both when
// base-calls match and when they mismatch. When they match, the two quality scores can contribute
// to a higher "consensus score". When they mismatch, the quality scores can settle which base-call
// should be retained and then contribute to a reduced quality score. In either case, information
// from the two base-calls, which is to say two independent data points, are put to use instead of
// being wasted.
mod correct;

// `errors` contains all the custom errors for the library.
pub mod errors;
pub use errors::{Error, Result};

// `workflow` contains a fluent interface using the "type state builder pattern" for library users
// to quickly get up and running with `libpairassembly`.
mod workflow;

// Symbols made available for library users are exposed in `prelude` and re-exported here.
pub mod prelude;

// re-exports
pub use prelude::*;
