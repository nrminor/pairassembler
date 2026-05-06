use std::{mem, time::Instant};

use color_eyre::eyre::{Result, WrapErr, bail};
use libpairassembly::{Assembler, OverlapSearch, OwnedSequenceRead, PairInput};
use noodles::fastq::Record as FastqRecord;
use rayon::prelude::*;
use tracing::info;

use crate::{
    RunRequest,
    progress::ProgressReporter,
    report::{self, RunContext, RunSummary},
    stats::{AssemblyStats, UnmergedReason},
};

use super::{
    IO_BATCH_SIZE,
    input::{PairedBatch, PairedFastqReaders},
    output::OutputHandles,
    records::{FastqMateKey, FastqReadView, mate_keys_match},
};

const MIN_PAIRS_PER_PARALLEL_TASK: usize = 128;

enum PairAssemblyOutcome {
    Merged(OwnedSequenceRead),
    Unmerged(UnmergedReason),
    SkippedMateMismatch,
}

pub(super) struct MergeOrchestrator<'request> {
    request: &'request RunRequest,
    inputs: PairedFastqReaders,
    outputs: OutputHandles,
    assembler: Assembler,
    no_correct: bool,
    no_validate: bool,
    context: RunContext,
    stats: AssemblyStats,
    progress: ProgressReporter,
    started_at: Instant,
    outcomes: Vec<Result<PairAssemblyOutcome>>,
}

impl<'request> MergeOrchestrator<'request> {
    pub(super) fn new(request: &'request RunRequest) -> Result<Self> {
        let context = RunContext::from_request(request);
        Ok(Self {
            request,
            inputs: PairedFastqReaders::spawn(
                &request.input1,
                &request.input2,
                context.input_label(),
            )?,
            outputs: OutputHandles::new(request)?,
            assembler: Assembler::builder()
                .with_overlap_params(request.settings.overlap_settings)
                .with_validator(request.settings.validation_settings)
                .build()?,
            no_correct: request.settings.no_correct,
            no_validate: request.settings.no_validate,
            context,
            stats: AssemblyStats::new(!request.settings.no_correct),
            progress: ProgressReporter::new(request.ui.progress_mode, request.progress_every),
            started_at: Instant::now(),
            outcomes: Vec::with_capacity(IO_BATCH_SIZE),
        })
    }

    pub(super) fn run(mut self) -> Result<()> {
        loop {
            let batch = match self.inputs.next_batch() {
                Ok(Some(batch)) => batch,
                Ok(None) => break,
                Err(error) => {
                    self.inputs.cancel();
                    return Err(error);
                },
            };

            let result = (|| -> Result<()> {
                self.record_active_batch(&batch)?;
                self.assemble_active_batch(&batch);
                self.write_outcomes(&batch)
            })();

            if let Err(error) = result {
                self.inputs.cancel();
                return Err(error);
            }

            self.inputs.recycle(batch);
        }

        self.finish()
    }

    fn assemble_active_batch(&mut self, batch: &PairedBatch) {
        let config = self.assembler.config().clone();
        let no_correct = self.no_correct;
        let no_validate = self.no_validate;
        let mut outcomes = mem::take(&mut self.outcomes);

        batch
            .r1_records()
            .par_iter()
            .zip(batch.r2_records().par_iter())
            .with_min_len(MIN_PAIRS_PER_PARALLEL_TASK)
            .map_init(
                || Assembler::from_config(config.clone()),
                |assembler, (r1, r2)| {
                    if !mate_keys_match(r1, r2) {
                        return Ok(PairAssemblyOutcome::SkippedMateMismatch);
                    }

                    let pair_id = FastqMateKey::from_record(r1).as_str()?;
                    let pair_input = PairInput::new(
                        FastqReadView::from_record(r1, pair_id)?,
                        FastqReadView::from_record(r2, pair_id)?,
                    );

                    let overlap = match assembler.on_pair(&pair_input)?.find_overlap()? {
                        OverlapSearch::Found(overlap) => overlap,
                        OverlapSearch::NoOverlap(_) => {
                            return Ok(PairAssemblyOutcome::Unmerged(
                                UnmergedReason::NoAcceptableOverlap,
                            ));
                        },
                    };

                    let read = if no_validate {
                        let merged = overlap.merge()?;
                        if no_correct {
                            merged.into_owned_read()?
                        } else {
                            merged.correct()?.into_owned_read()?
                        }
                    } else {
                        let Ok(validated) = overlap.validate() else {
                            return Ok(PairAssemblyOutcome::Unmerged(
                                UnmergedReason::OverlapRejectedByValidation,
                            ));
                        };

                        let merged = validated.merge()?;
                        if no_correct {
                            merged.into_owned_read()?
                        } else {
                            merged.correct()?.into_owned_read()?
                        }
                    };

                    Ok(PairAssemblyOutcome::Merged(read))
                },
            )
            .collect_into_vec(&mut outcomes);

        self.outcomes = outcomes;
    }

    fn record_pair_seen(&mut self, r1: &FastqRecord, r2: &FastqRecord) {
        self.stats
            .record_pair_seen(r1.sequence().len(), r2.sequence().len());
    }

    fn record_active_batch(&mut self, batch: &PairedBatch) -> Result<()> {
        for (r1, r2) in batch.records() {
            self.record_pair_seen(r1, r2);

            if self.record_mate_mismatch_if_needed(r1, r2)? {
                self.progress.maybe_report(&self.stats);
            }
        }

        Ok(())
    }

    fn record_mate_mismatch_if_needed(
        &mut self,
        r1: &FastqRecord,
        r2: &FastqRecord,
    ) -> Result<bool> {
        if mate_keys_match(r1, r2) {
            return Ok(false);
        }

        self.stats.record_mate_id_mismatch();
        if self.stats.mate_id_mismatches > self.request.settings.max_mate_id_mismatches {
            bail!(
                "paired FASTQ inputs appear to be in different orders\nsource: {}\nmate_id_mismatches: {}\nmax_mate_id_mismatches: {}\nlast_r1_header: {}\nlast_r2_header: {}\npairs_seen_before_failure: {}\nhelp: pairasm expects R1 and R2 FASTQs to be sorted in the same pair order; repair or re-sort pairing before running pairasm",
                self.context.input_label(),
                self.stats.mate_id_mismatches,
                self.request.settings.max_mate_id_mismatches,
                String::from_utf8_lossy(r1.name().as_ref()),
                String::from_utf8_lossy(r2.name().as_ref()),
                self.stats.pairs_seen,
            );
        }

        Ok(true)
    }

    fn write_outcomes(&mut self, batch: &PairedBatch) -> Result<()> {
        for ((r1, r2), outcome) in batch.records().zip(self.outcomes.drain(..)) {
            match outcome? {
                PairAssemblyOutcome::Merged(read) => {
                    self.outputs.write_merged(&read)?;
                    self.stats.record_merged(read.sequence_bytes().len());
                },
                PairAssemblyOutcome::Unmerged(reason) => {
                    let was_written = self.outputs.write_unmerged_records(r1, r2)?;
                    self.stats.record_unmerged(reason, was_written);
                },
                PairAssemblyOutcome::SkippedMateMismatch => {},
            }

            self.progress.maybe_report(&self.stats);
        }

        Ok(())
    }

    fn finish(self) -> Result<()> {
        self.progress.finish();

        let output_result = self.outputs.finish();
        let input_result = self.inputs.finish();
        output_result?;
        input_result?;

        let summary = RunSummary::from_stats(self.context, self.stats, self.started_at.elapsed());
        if self.request.ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = self.request.summary.as_deref() {
            report::write_summary_json(path, &summary).wrap_err_with(|| {
                format!(
                    "failed to write JSON run summary: {}\nhelp: check that the parent directory exists and is writable",
                    path.display()
                )
            })?;
        }

        info!(
            pairs_seen = summary.stats.pairs_seen,
            pairs_processed = summary.stats.pairs_processed,
            pairs_merged = summary.stats.pairs_merged,
            pairs_unmerged = summary.stats.pairs_unmerged,
            mate_id_mismatches = summary.stats.mate_id_mismatches,
            elapsed_seconds = summary.elapsed_seconds,
            outcome = "completed",
            "pairasm_run_finished"
        );

        Ok(())
    }
}
