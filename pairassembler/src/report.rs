use std::{fmt::Write as _, fs::File, io::BufWriter, path::Path, time::Duration};

use color_eyre::eyre::Result;
use serde::Serialize;

use crate::RunRequest;
use crate::stats::AssemblyStats;

/// Input and configuration context attached to a run summary.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RunContext {
    pub input1: String,
    pub input2: String,
    pub output_file: Option<String>,
    pub unmerged_out: Option<String>,
    pub max_mate_id_mismatches: u64,
}

impl RunContext {
    /// Build run context from the resolved CLI request.
    #[must_use]
    pub fn from_request(request: &RunRequest) -> Self {
        Self {
            input1: request.input1.clone(),
            input2: request.input2.clone(),
            output_file: request.output_file.clone(),
            unmerged_out: request.unmerged_output.clone(),
            max_mate_id_mismatches: request.settings.max_mate_id_mismatches,
        }
    }

    /// Return a compact label for diagnostic messages.
    #[must_use]
    pub fn input_label(&self) -> String {
        format!("local-paired:{}|{}", self.input1, self.input2)
    }
}

/// Compact owned summary for stderr presentation and JSON serialization.
#[derive(Clone, Debug, Serialize)]
pub struct RunSummary {
    pub context: RunContext,
    pub elapsed_seconds: f64,
    pub pairs_per_second: f64,
    pub merge_fraction: f64,
    pub unmerged_fraction: f64,
    pub stats: AssemblyStats,
}

impl RunSummary {
    /// Build an owned summary from run context, counters, and elapsed wall time.
    #[must_use]
    pub fn from_stats(context: RunContext, stats: AssemblyStats, elapsed: Duration) -> Self {
        let elapsed_seconds = elapsed.as_secs_f64();
        Self {
            context,
            pairs_per_second: rate(stats.pairs_seen, elapsed_seconds),
            merge_fraction: fraction(stats.pairs_merged, stats.pairs_processed),
            unmerged_fraction: fraction(stats.pairs_unmerged, stats.pairs_processed),
            elapsed_seconds,
            stats,
        }
    }
}

/// Print a human-readable run summary to stderr.
pub fn print_summary(summary: &RunSummary) {
    eprint!("{}", render_summary(summary));
}

/// Write a JSON summary to disk.
///
/// # Errors
///
/// Returns an error when the file cannot be created or the JSON cannot be serialized.
pub fn write_summary_json(path: &Path, summary: &RunSummary) -> Result<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, summary)?;
    Ok(())
}

fn render_summary(summary: &RunSummary) -> String {
    let mut output = String::new();
    output.push('\n');
    output.push_str("pairassembler summary\n");
    write_context(&mut output, summary);
    write_totals(&mut output, summary);
    output
}

fn write_context(output: &mut String, summary: &RunSummary) {
    let _ = writeln!(output, "  input 1:             {}", summary.context.input1);
    let _ = writeln!(output, "  input 2:             {}", summary.context.input2);
    if let Some(path) = &summary.context.output_file {
        let _ = writeln!(output, "  merged output:       {path}");
    }
    if let Some(path) = &summary.context.unmerged_out {
        let _ = writeln!(output, "  unmerged output:     {path}");
    }
}

fn write_totals(output: &mut String, summary: &RunSummary) {
    let stats = &summary.stats;
    let _ = writeln!(
        output,
        "  elapsed:             {:.2}s",
        summary.elapsed_seconds
    );
    let _ = writeln!(
        output,
        "  throughput:          {:.1} pairs/s",
        summary.pairs_per_second
    );
    let _ = writeln!(output, "  pairs seen:          {}", stats.pairs_seen);
    let _ = writeln!(output, "  pairs processed:     {}", stats.pairs_processed);
    let _ = writeln!(
        output,
        "  pairs merged:        {} ({:.2}%)",
        stats.pairs_merged,
        summary.merge_fraction * 100.0
    );
    let _ = writeln!(
        output,
        "  pairs unmerged:      {} ({:.2}%)",
        stats.pairs_unmerged,
        summary.unmerged_fraction * 100.0
    );
    let _ = writeln!(
        output,
        "    no overlap:        {}",
        stats.pairs_unmerged_no_overlap
    );
    let _ = writeln!(
        output,
        "    validation reject: {}",
        stats.pairs_unmerged_validation_rejected
    );
    let _ = writeln!(
        output,
        "  unmerged written:    {}",
        stats.unmerged_pairs_written
    );
    let _ = writeln!(
        output,
        "  mate ID mismatches:  {}",
        stats.mate_id_mismatches
    );
    let _ = writeln!(output, "  bases in:            {}", stats.bases_in);
    let _ = writeln!(output, "  bases merged:        {}", stats.bases_merged);
    let _ = writeln!(
        output,
        "  quality correction:  {}",
        if stats.quality_correction_enabled {
            "enabled"
        } else {
            "disabled"
        }
    );
}

fn fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        u64_to_f64(numerator) / u64_to_f64(denominator)
    }
}

fn rate(total: u64, elapsed_seconds: f64) -> f64 {
    if elapsed_seconds <= f64::EPSILON {
        0.0
    } else {
        u64_to_f64(total) / elapsed_seconds
    }
}

fn u64_to_f64(value: u64) -> f64 {
    value.to_string().parse::<f64>().unwrap_or(f64::INFINITY)
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::{RunContext, RunSummary, render_summary};
    use crate::stats::{AssemblyStats, UnmergedReason};

    #[test]
    fn render_summary_describes_unmerged_as_normal_outcome() {
        let mut stats = AssemblyStats::new(true);
        stats.record_pair_seen(4, 4);
        stats.record_unmerged(UnmergedReason::NoAcceptableOverlap, false);
        let summary = RunSummary::from_stats(
            RunContext {
                input1: "r1.fastq".to_owned(),
                input2: "r2.fastq".to_owned(),
                output_file: None,
                unmerged_out: None,
                max_mate_id_mismatches: 3,
            },
            stats,
            Duration::from_secs(1),
        );

        let rendered = render_summary(&summary);
        assert!(rendered.contains("pairs unmerged"));
        assert!(rendered.contains("no overlap"));
    }
}
