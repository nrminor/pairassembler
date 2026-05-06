use std::path::Path;

use color_eyre::eyre::Result;

use crate::RunRequest;

mod input;
mod orchestrator;
mod output;
mod records;

const IO_BATCH_SIZE: usize = 8192;

// Two reusable batches per mate are enough to let reader threads fill the next
// batch while the current batch is being assembled and written. Larger pools did
// not improve the default benchmarks enough to justify the extra resident memory.
const READ_BATCH_POOL_SIZE: usize = 2;

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
