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
use libpairassembly::{OverlapParams, OverlapValidator};
use pairassembler::{RunRequest, RunSettings, cli::Cli, merging};

fn main() -> Result<()> {
    color_eyre::install()?;
    let cli = Cli::parse();
    cli.init_tracing()?;
    tracing::info!(version = env!("CARGO_PKG_VERSION"), "starting pairasm");
    let ui = cli.ui_policy();

    let Cli {
        verbosity: _,
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
        no_correct,
        progress_every,
        summary,
        max_mate_id_mismatches,
    } = cli;

    let overlap_settings = OverlapParams::default()
        .with_overlap_diff_max(overlap_diff_max)
        .with_min_overlap(min_overlap)
        .with_diff_percent_max(diff_percent_max)
        .with_min_comparisons(min_comparisons);
    let validation_settings = OverlapValidator::default()
        .with_k(k)
        .with_min_complexity_score(min_complexity_score);
    let settings = RunSettings::new(
        overlap_settings,
        validation_settings,
        no_correct,
        max_mate_id_mismatches,
    );
    let request = RunRequest {
        input1,
        input2,
        output_file,
        unmerged_output: unmerged_out,
        summary,
        progress_every,
        ui,
        settings,
    };
    merging::run(&request)?;

    Ok(())
}
