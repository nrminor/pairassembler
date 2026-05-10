mod artifacts;
mod cli;
mod commands;
mod config;
mod db;
mod fastq;
mod fetch;
mod model;
mod prepare;
mod process;
mod products;
mod report;
mod run;
mod shell;
mod ui;
mod validate;
mod vcs;

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
        BenchCommand::Report(options) => report::report(&options),
        BenchCommand::WorkflowPhase(options) => {
            ui::print_workflow_phase(&options.step, &options.title);
            Ok(())
        },
    }
}
