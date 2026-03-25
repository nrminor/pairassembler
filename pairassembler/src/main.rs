// dev allowances
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]
//
// crate-level lints
#![warn(
    clippy::pedantic,
    clippy::perf,
    // clippy::todo,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::complexity,
    clippy::correctness,
    clippy::absolute_paths,
    clippy::style
)]

use clap::Parser;
use color_eyre::{self, Result};
use pairassembler::{
    RunSettings,
    cli::{
        self, Cli,
        Commands::{Correct, Merge},
    },
    merging,
};
use tracing_subscriber::EnvFilter;

#[tokio::main]
async fn main() -> Result<()> {
    // set up the color-eyre display and tracer
    utils::setup()?;

    // parse command line arguments and handle verbosity and strictness, which will essentially be globally
    // scoped settings
    let Cli {
        verbose,
        command,
        strict,
    } = Cli::parse();

    // match on a provided command if given, printing info if not
    match command {
        // Run paired read merging with user settings
        Some(Merge {
            input1,
            input2,
            output_file,
            unmerged_out,
            no_correct,
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            k,
            min_entropy,
        }) => {
            let settings = RunSettings::new(
                overlap_diff_max,
                min_overlap,
                diff_percent_max,
                min_comparisons,
                k,
                min_entropy,
                no_correct,
            );
            merging::run(input1, input2, output_file, unmerged_out, settings).await?;
        },

        // Run correction but not merging with user settings
        Some(Correct {
            input1,
            input2,
            output_file,
            unmerged_out,
            overlap_diff_max,
            min_overlap,
            diff_percent_max,
            min_comparisons,
            k,
            min_entropy,
        }) => {
            todo!()
        },

        // Run validation that reads are in an appropriate in one or two provided input files
        Some(cli::Commands::Validate { input1, input2 }) => {
            todo!()
        },

        // Sort input files by read ID so that they can be run through the above command successfully
        Some(cli::Commands::Sort { input1, input2 }) => {
            todo!()
        },

        // No subcommand provided and thus no computation is requested. Closing with the application's info.
        None => {
            eprintln!("{}\n", cli::INFO);
        },
    }

    Ok(())
}

mod utils {
    use std::env;

    use tracing_subscriber::fmt;

    use super::*;

    pub(super) fn setup() -> Result<()> {
        if env::var("RUST_LIB_BACKTRACE").is_err() {
            // UNSAFE: temporarily allowing until we find a better solution; this shouldn't
            // be on in a released library anyway
            unsafe { env::set_var("RUST_LIB_BACKTRACE", "1") }
        }
        color_eyre::install()?;

        if env::var("RUST_LOG").is_err() {
            // UNSAFE: temporarily allowing until we find a better solution; this shouldn't
            // be on in a released library anyway
            unsafe { env::set_var("RUST_LOG", "info") }
        }
        fmt::fmt()
            .with_env_filter(EnvFilter::from_default_env())
            .init();

        Ok(())
    }
}
