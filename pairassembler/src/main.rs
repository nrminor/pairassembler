#![allow(dead_code, unused_imports, unused_variables, unused_mut)]
#![warn(
    clippy::pedantic,
    clippy::perf,
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
    utils::setup()?;

    let Cli {
        verbose,
        command,
        strict,
    } = Cli::parse();

    match command {
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
            min_complexity_score,
        }) => {
            let settings = RunSettings::new(
                overlap_diff_max,
                min_overlap,
                diff_percent_max,
                min_comparisons,
                k,
                min_complexity_score,
                no_correct,
            );
            merging::run(input1, input2, output_file, unmerged_out, settings).await?;
        },

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
            min_complexity_score,
        }) => {
            todo!()
        },

        Some(cli::Commands::Validate { input1, input2 }) => {
            todo!()
        },

        Some(cli::Commands::Sort { input1, input2 }) => {
            todo!()
        },

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
            // SAFETY: process environment defaults are set during single-threaded startup before
            // worker tasks are spawned.
            unsafe { env::set_var("RUST_LIB_BACKTRACE", "1") }
        }
        color_eyre::install()?;

        if env::var("RUST_LOG").is_err() {
            // SAFETY: process environment defaults are set during single-threaded startup before
            // worker tasks are spawned.
            unsafe { env::set_var("RUST_LOG", "info") }
        }
        fmt::fmt()
            .with_env_filter(EnvFilter::from_default_env())
            .init();

        Ok(())
    }
}
