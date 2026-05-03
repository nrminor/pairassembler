use std::time::Instant;

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

use crate::stats::AssemblyStats;

/// Progress rendering mode selected from quietness and terminal capabilities.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressMode {
    /// Continuously refresh a single stderr status line.
    Live,
    /// Emit plain periodic stderr lines without terminal cursor control.
    Plain,
    /// Disable progress reporting entirely.
    Off,
}

/// Emits progress updates after a configurable number of complete pairs have been read.
pub struct ProgressReporter {
    mode: ProgressMode,
    every: u64,
    next_threshold: u64,
    started_at: Instant,
    bar: Option<ProgressBar>,
}

impl ProgressReporter {
    /// Construct a reporter that updates every `every` seen pairs in the selected mode.
    #[must_use]
    pub fn new(mode: ProgressMode, every: u64) -> Self {
        let bar = (mode == ProgressMode::Live).then(|| {
            let bar = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
            bar.set_style(
                ProgressStyle::default_spinner()
                    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                    .template("{spinner} {msg}")
                    .unwrap_or_else(|_| ProgressStyle::default_spinner()),
            );
            bar
        });

        Self {
            mode,
            every,
            next_threshold: every,
            started_at: Instant::now(),
            bar,
        }
    }

    /// Emit a progress update if the next threshold has been crossed.
    pub fn maybe_report(&mut self, stats: &AssemblyStats) {
        if self.mode == ProgressMode::Off
            || self.every == 0
            || stats.pairs_seen < self.next_threshold
        {
            return;
        }

        let message = render_message(stats, self.started_at.elapsed().as_secs_f64());

        match self.mode {
            ProgressMode::Live => {
                if let Some(bar) = &self.bar {
                    bar.set_message(message);
                    bar.tick();
                }
            },
            ProgressMode::Plain => eprintln!("{message}"),
            ProgressMode::Off => {},
        }

        self.next_threshold = self.next_threshold.saturating_add(self.every);
    }

    /// Clear any live terminal progress rendering before final summary output.
    pub fn finish(&self) {
        if let Some(bar) = &self.bar {
            bar.finish_and_clear();
        }
    }
}

fn render_message(stats: &AssemblyStats, elapsed_seconds: f64) -> String {
    let top_unmerged = stats.top_unmerged_reason().map_or_else(
        || "none".to_owned(),
        |(reason, count)| format!("{reason}:{count}"),
    );

    format!(
        "pairs={} merged={} unmerged={} ({:.1}% merged) mate_mismatches={} {:.0} pairs/s top_unmerged={}",
        stats.pairs_seen,
        stats.pairs_merged,
        stats.pairs_unmerged,
        fraction(stats.pairs_merged, stats.pairs_processed) * 100.0,
        stats.mate_id_mismatches,
        rate(stats.pairs_seen, elapsed_seconds),
        top_unmerged,
    )
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
    use super::{ProgressMode, ProgressReporter, rate};

    #[test]
    fn rate_is_zero_when_elapsed_is_zero() {
        assert!((rate(42, 0.0) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn live_reporter_constructs_progress_bar() {
        let reporter = ProgressReporter::new(ProgressMode::Live, 100);
        assert_eq!(reporter.mode, ProgressMode::Live);
        assert!(reporter.bar.is_some());
    }
}
