mod cli;
mod commands;
mod config;
mod fastq;
mod fetch;
mod model;
mod prepare;
mod process;
mod run;
mod shell;
mod summarize;
mod validate;

use clap::Parser;
use color_eyre::eyre::Result;

use crate::cli::{BenchCommand, Cli};

fn main() -> Result<()> {
    color_eyre::install()?;
    let cli = Cli::parse();

    match cli.command {
        BenchCommand::Check => config::check_tools(),
        BenchCommand::Fetch(options) => fetch::fetch_ena(&options),
        BenchCommand::Prepare(options) => prepare::prepare_subsets(&options),
        BenchCommand::Run(options) => run::run_matrix(&options),
        BenchCommand::Summarize(options) => summarize::summarize(&options),
    }
}
