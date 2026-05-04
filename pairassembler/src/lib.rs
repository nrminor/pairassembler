#![warn(
    clippy::perf,
    clippy::unwrap_used,
    clippy::complexity,
    clippy::correctness,
    clippy::absolute_paths,
    clippy::style
)]

use std::path::PathBuf;

use libpairassembly::{OverlapParams, OverlapValidator};

use crate::cli::UiPolicy;

pub mod cli;
pub mod merging;
pub mod progress;
pub mod report;
pub mod stats;

#[derive(Debug)]
pub struct RunSettings {
    no_correct: bool,
    max_mate_id_mismatches: u64,
    overlap_settings: OverlapParams,
    validation_settings: OverlapValidator,
}

impl RunSettings {
    #[must_use]
    pub const fn new(
        overlap_settings: OverlapParams,
        validation_settings: OverlapValidator,
        no_correct: bool,
        max_mate_id_mismatches: u64,
    ) -> Self {
        RunSettings {
            no_correct,
            max_mate_id_mismatches,
            overlap_settings,
            validation_settings,
        }
    }
}

#[derive(Debug)]
pub struct RunRequest {
    pub input1: String,
    pub input2: String,
    pub output_file: Option<String>,
    pub unmerged_output: Option<String>,
    pub summary: Option<PathBuf>,
    pub progress_every: u64,
    pub ui: UiPolicy,
    pub settings: RunSettings,
}
