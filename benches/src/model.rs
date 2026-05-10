use std::{num::NonZeroUsize, path::PathBuf};

use clap::ValueEnum;
use serde::Deserialize;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum Tool {
    Pairasm,
    Fastp,
    Bbmerge,
    Vsearch,
}

impl Tool {
    fn as_str(self) -> &'static str {
        match self {
            Tool::Pairasm => "pairasm",
            Tool::Fastp => "fastp",
            Tool::Bbmerge => "bbmerge",
            Tool::Vsearch => "vsearch",
        }
    }
}

impl std::fmt::Display for Tool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Clone, Debug)]
pub struct Dataset {
    pub name: String,
    pub accession: String,
    pub read_pair_cap: Option<NonZeroUsize>,
    pub note: String,
}

#[derive(Clone, Debug)]
pub struct ToolPaths {
    pub pairasm: PathBuf,
    pub fastp: PathBuf,
    pub bbmerge: PathBuf,
    pub vsearch: PathBuf,
    pub hyperfine: PathBuf,
}

#[derive(Debug)]
pub struct SourceMetadata {
    pub name: String,
    pub accession: String,
    pub r1: PathBuf,
    pub r2: PathBuf,
}

#[derive(Debug)]
pub struct SubsetMetadata {
    pub name: String,
    pub accession: String,
    pub read_pairs: usize,
    pub r1: PathBuf,
    pub r2: PathBuf,
}

#[derive(Debug)]
pub struct ToolCommand {
    pub tool: Tool,
    pub args: Vec<String>,
    pub merged_output: PathBuf,
}

#[derive(Debug, Deserialize)]
pub struct HyperfineReport {
    pub results: Vec<HyperfineResult>,
}

#[derive(Debug, Deserialize)]
pub struct HyperfineResult {
    pub mean: f64,
    pub stddev: Option<f64>,
    pub median: f64,
    pub min: f64,
    pub max: f64,
    pub user: f64,
    pub system: f64,
}
