use std::path::Path;

use color_eyre::eyre::Result;

use crate::RunRequest;

mod input;
mod orchestrator;
mod output;
mod records;

const IO_BATCH_SIZE: usize = 8192;

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|extension| extension == "gz")
}

/// Run pair merging over two FASTQ inputs.
///
/// # Errors
///
/// Returns an error when input files cannot be read, paired inputs violate the ordering
/// contract too often, output cannot be written, or a non-biological assembly invariant fails.
pub fn run(request: &RunRequest) -> Result<()> {
    orchestrator::MergeOrchestrator::new(request)?.run()
}
