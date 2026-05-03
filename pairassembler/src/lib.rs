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

pub mod merging {
    use std::{
        fs::File,
        io::{self, BufRead, BufReader, BufWriter, Write},
        path::Path,
        str,
        time::Instant,
    };

    use color_eyre::eyre::{Result, WrapErr, bail};
    use flate2::{Compression, read::GzDecoder, write::GzEncoder};
    use libpairassembly::{Assembler, OverlapSearch, OwnedSequenceRead, PairInput, SequenceRead};
    use noodles::fastq::{
        Reader as FastqReader, Record as FastqRecord, io::Writer as FastqWriter, record::Definition,
    };
    use tokio::task;
    use tracing::info;

    use crate::{
        RunRequest, RunSettings,
        progress::ProgressReporter,
        report::{self, RunContext, RunSummary},
        stats::{AssemblyStats, UnmergedReason},
    };

    enum PairAssemblyOutcome {
        Merged(OwnedSequenceRead),
        Unmerged(UnmergedReason),
    }

    enum FastqOutput {
        Plain(FastqWriter<Box<dyn Write + Send>>),
        Gzip(FastqWriter<GzEncoder<Box<dyn Write + Send>>>),
    }

    impl FastqOutput {
        fn new(path: Option<&str>) -> Result<Self> {
            match path {
                Some(path) if is_gzip_path(Path::new(path)) => {
                    let file = File::create(path).wrap_err_with(|| {
                        format!(
                            "failed to create gzip FASTQ output: {path}\nhelp: check that the parent directory exists and is writable"
                        )
                    })?;
                    let writer: Box<dyn Write + Send> = Box::new(BufWriter::new(file));
                    Ok(Self::Gzip(FastqWriter::new(GzEncoder::new(
                        writer,
                        Compression::default(),
                    ))))
                },
                Some(path) => {
                    let file = File::create(path).wrap_err_with(|| {
                        format!(
                            "failed to create FASTQ output: {path}\nhelp: check that the parent directory exists and is writable"
                        )
                    })?;
                    Ok(Self::Plain(FastqWriter::new(Box::new(BufWriter::new(
                        file,
                    )))))
                },
                None => Ok(Self::Plain(FastqWriter::new(Box::new(BufWriter::new(
                    io::stdout(),
                ))))),
            }
        }

        fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
            match self {
                Self::Plain(writer) => writer.write_record(record)?,
                Self::Gzip(writer) => writer.write_record(record)?,
            }
            Ok(())
        }

        fn finish(self) -> Result<()> {
            match self {
                Self::Plain(mut writer) => writer.get_mut().flush()?,
                Self::Gzip(writer) => {
                    let encoder = writer.into_inner();
                    let mut inner = encoder.finish()?;
                    inner.flush()?;
                },
            }
            Ok(())
        }
    }

    struct OutputHandles {
        merged: FastqOutput,
        unmerged: Option<FastqOutput>,
    }

    impl OutputHandles {
        fn new(request: &RunRequest) -> Result<Self> {
            Ok(Self {
                merged: FastqOutput::new(request.output_file.as_deref())?,
                unmerged: match request.unmerged_output.as_deref() {
                    Some(path) => Some(FastqOutput::new(Some(path))?),
                    None => None,
                },
            })
        }

        fn finish(self) -> Result<()> {
            self.merged
                .finish()
                .wrap_err("failed to finalize merged FASTQ output")?;
            if let Some(output) = self.unmerged {
                output
                    .finish()
                    .wrap_err("failed to finalize unmerged FASTQ output")?;
            }
            Ok(())
        }
    }

    fn open_fastq_reader(path: impl AsRef<Path>) -> Result<FastqReader<Box<dyn BufRead>>> {
        let path = path.as_ref();
        let file_handle = File::open(path).wrap_err_with(|| {
            format!(
                "failed to open FASTQ input: {}\nhelp: check that the input path exists and is readable",
                path.display()
            )
        })?;
        let read_buffer = BufReader::new(file_handle);

        let reader: Box<dyn BufRead> = if is_gzip_path(path) {
            Box::new(BufReader::new(GzDecoder::new(read_buffer)))
        } else {
            Box::new(read_buffer)
        };
        Ok(FastqReader::new(reader))
    }

    fn is_gzip_path(path: &Path) -> bool {
        path.extension().is_some_and(|extension| extension == "gz")
    }

    fn sequence_read_from_record<'a>(
        record: &'a FastqRecord,
        pair_id: &'a str,
    ) -> Result<SequenceRead<'a>> {
        let seq =
            str::from_utf8(record.sequence()).wrap_err("FASTQ sequence is not valid UTF-8")?;
        let qual = str::from_utf8(record.quality_scores())
            .wrap_err("FASTQ quality string is not valid UTF-8")?;

        Ok(SequenceRead::try_new(pair_id, seq, qual)?)
    }

    fn mate_key(header: &[u8]) -> &[u8] {
        let first_token = header
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(header);

        match first_token {
            [prefix @ .., b'/', b'1' | b'2'] => prefix,
            _ => first_token,
        }
    }

    fn merged_record(read: &OwnedSequenceRead) -> FastqRecord {
        FastqRecord::new(
            Definition::new(read.id(), ""),
            read.sequence_bytes(),
            read.quality_bytes(),
        )
    }

    fn process_pair(
        assembler: &Assembler,
        fwd: &FastqRecord,
        rev: &FastqRecord,
        settings: &RunSettings,
    ) -> Result<PairAssemblyOutcome> {
        let pair_key = str::from_utf8(mate_key(fwd.name().as_ref()))
            .wrap_err("FASTQ read name is not valid UTF-8")?;
        let read1 = sequence_read_from_record(fwd, pair_key)?;
        let read2 = sequence_read_from_record(rev, pair_key)?;
        let pair_input = PairInput::new(read1, read2);

        let overlap = match assembler.on_pair(&pair_input)?.find_overlap()? {
            OverlapSearch::Found(ctx) => ctx,
            OverlapSearch::NoOverlap(_) => {
                return Ok(PairAssemblyOutcome::Unmerged(
                    UnmergedReason::NoAcceptableOverlap,
                ));
            },
        };

        let Ok(validated) = overlap.validate() else {
            return Ok(PairAssemblyOutcome::Unmerged(
                UnmergedReason::OverlapRejectedByValidation,
            ));
        };
        let merged = validated.merge()?;
        let read = if settings.no_correct {
            merged.into_owned_read()?
        } else {
            merged.correct()?.into_owned_read()?
        };

        Ok(PairAssemblyOutcome::Merged(read))
    }

    /// Run pair merging over two FASTQ inputs.
    ///
    /// # Errors
    ///
    /// Returns an error when input files cannot be read, paired inputs violate the ordering
    /// contract too often, output cannot be written, or a non-biological assembly invariant fails.
    pub async fn run(request: RunRequest) -> Result<()> {
        task::spawn_blocking(move || run_sync(&request))
            .await
            .wrap_err("merge worker task failed to join")?
    }

    fn run_sync(request: &RunRequest) -> Result<()> {
        let mut fastq_reader1 = open_fastq_reader(&request.input1)?;
        let mut fastq_reader2 = open_fastq_reader(&request.input2)?;
        let mut records1 = fastq_reader1.records();
        let mut records2 = fastq_reader2.records();
        let mut outputs = OutputHandles::new(request)?;
        let assembler = Assembler::builder()
            .with_overlap_params(request.settings.overlap_settings)
            .with_validator(request.settings.validation_settings)
            .build()?;
        let context = RunContext::from_request(request);
        let mut stats = AssemblyStats::new(!request.settings.no_correct);
        let mut progress = ProgressReporter::new(request.ui.progress_mode, request.progress_every);
        let started_at = Instant::now();

        loop {
            let Some((fwd, rev)) =
                next_pair(records1.next(), records2.next(), request, &context, &stats)?
            else {
                break;
            };

            stats.record_pair_seen(fwd.sequence().len(), rev.sequence().len());

            if handle_mate_mismatch(&fwd, &rev, request, &context, &mut stats)? {
                progress.maybe_report(&stats);
                continue;
            }

            let outcome = process_pair(&assembler, &fwd, &rev, &request.settings)?;
            write_outcome(outcome, &fwd, &rev, &mut outputs, &mut stats)?;

            progress.maybe_report(&stats);
        }

        progress.finish();
        outputs.finish()?;

        let summary = RunSummary::from_stats(context, stats, started_at.elapsed());
        if request.ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = request.summary.as_deref() {
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

    fn next_pair(
        next1: Option<io::Result<FastqRecord>>,
        next2: Option<io::Result<FastqRecord>>,
        request: &RunRequest,
        context: &RunContext,
        stats: &AssemblyStats,
    ) -> Result<Option<(FastqRecord, FastqRecord)>> {
        match (next1, next2) {
            (Some(Ok(fwd)), Some(Ok(rev))) => Ok(Some((fwd, rev))),
            (None, None) => Ok(None),
            (Some(Err(error)), _) => {
                Err(error).wrap_err_with(|| format!("failed to read {}", request.input1))
            },
            (_, Some(Err(error))) => {
                Err(error).wrap_err_with(|| format!("failed to read {}", request.input2))
            },
            (Some(Ok(_)), None) | (None, Some(Ok(_))) => {
                bail!(
                    "paired FASTQ inputs have different record counts\nsource: {}\ncomplete_pairs_seen: {}\nhelp: pairasm expects R1 and R2 FASTQs to be in the same order and have the same number of records",
                    context.input_label(),
                    stats.pairs_seen,
                );
            },
        }
    }

    fn handle_mate_mismatch(
        fwd: &FastqRecord,
        rev: &FastqRecord,
        request: &RunRequest,
        context: &RunContext,
        stats: &mut AssemblyStats,
    ) -> Result<bool> {
        if mate_key(fwd.name().as_ref()) == mate_key(rev.name().as_ref()) {
            return Ok(false);
        }

        stats.record_mate_id_mismatch();
        if stats.mate_id_mismatches > request.settings.max_mate_id_mismatches {
            bail!(
                "paired FASTQ inputs appear to be in different orders\nsource: {}\nmate_id_mismatches: {}\nmax_mate_id_mismatches: {}\nlast_r1_header: {}\nlast_r2_header: {}\npairs_seen_before_failure: {}\nhelp: pairasm expects R1 and R2 FASTQs to be sorted in the same pair order; repair or re-sort pairing before running pairasm",
                context.input_label(),
                stats.mate_id_mismatches,
                request.settings.max_mate_id_mismatches,
                String::from_utf8_lossy(fwd.name().as_ref()),
                String::from_utf8_lossy(rev.name().as_ref()),
                stats.pairs_seen,
            );
        }

        Ok(true)
    }

    fn write_outcome(
        outcome: PairAssemblyOutcome,
        fwd: &FastqRecord,
        rev: &FastqRecord,
        outputs: &mut OutputHandles,
        stats: &mut AssemblyStats,
    ) -> Result<()> {
        match outcome {
            PairAssemblyOutcome::Merged(read) => {
                stats.record_merged(read.sequence_bytes().len());
                outputs
                    .merged
                    .write_record(&merged_record(&read))
                    .wrap_err_with(|| {
                        format!(
                            "failed to write merged FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                            read.id()
                        )
                    })?;
            },
            PairAssemblyOutcome::Unmerged(reason) => {
                let was_written = if let Some(output) = &mut outputs.unmerged {
                    output.write_record(fwd)?;
                    output.write_record(rev)?;
                    true
                } else {
                    false
                };
                stats.record_unmerged(reason, was_written);
            },
        }

        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use libpairassembly::{OverlapParams, OverlapValidator};

        use super::{Definition, FastqRecord, handle_mate_mismatch, mate_key};
        use crate::{
            RunRequest, RunSettings, cli::UiPolicy, progress::ProgressMode, report::RunContext,
            stats::AssemblyStats,
        };

        #[test]
        fn mate_key_strips_slash_mate_suffix() {
            assert_eq!(mate_key(b"read123/1"), b"read123");
            assert_eq!(mate_key(b"read123/2"), b"read123");
            assert_eq!(mate_key(b"read123"), b"read123");
        }

        #[test]
        fn mate_key_uses_first_whitespace_token() {
            assert_eq!(mate_key(b"read123/1 instrument stuff"), b"read123");
            assert_eq!(mate_key(b"read123 comment"), b"read123");
        }

        #[test]
        fn mate_mismatch_threshold_fails_fast_without_counting_unmerged() {
            let request = RunRequest {
                input1: "r1.fastq".to_owned(),
                input2: "r2.fastq".to_owned(),
                output_file: None,
                unmerged_output: None,
                summary: None,
                progress_every: 100_000,
                ui: UiPolicy {
                    log_level: None,
                    show_summary: false,
                    progress_mode: ProgressMode::Off,
                },
                settings: RunSettings::new(
                    OverlapParams::default(),
                    OverlapValidator::default(),
                    false,
                    0,
                ),
            };
            let context = RunContext::from_request(&request);
            let fwd = FastqRecord::new(Definition::new("read-a/1", ""), "AAAA", "IIII");
            let rev = FastqRecord::new(Definition::new("read-b/2", ""), "TTTT", "IIII");
            let mut stats = AssemblyStats::new(true);
            stats.record_pair_seen(fwd.sequence().len(), rev.sequence().len());

            let result = handle_mate_mismatch(&fwd, &rev, &request, &context, &mut stats);

            assert!(result.is_err());
            assert_eq!(stats.mate_id_mismatches, 1);
            assert_eq!(stats.pairs_unmerged, 0);
            assert_eq!(stats.pairs_processed, 0);
        }
    }
}
