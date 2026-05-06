use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    mem,
    path::Path,
    str,
    sync::mpsc::{self, Receiver, SyncSender},
    thread::{self, JoinHandle},
    time::Instant,
};

use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use libpairassembly::{Assembler, OverlapSearch, OwnedSequenceRead, PairInput, SeqRecordView};
use noodles::fastq::{
    Reader as FastqReader, Record as FastqRecord, io::Writer as FastqWriter, record::Definition,
};
use rayon::prelude::*;
use tracing::info;

use crate::{
    RunRequest,
    progress::ProgressReporter,
    report::{self, RunContext, RunSummary},
    stats::{AssemblyStats, UnmergedReason},
};

const IO_BATCH_SIZE: usize = 8192;
const MIN_PAIRS_PER_PARALLEL_TASK: usize = 128;

enum PairAssemblyOutcome {
    Merged(OwnedSequenceRead),
    Unmerged(UnmergedReason),
    SkippedMateMismatch,
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

#[derive(Clone, Copy)]
enum Mate {
    R1,
    R2,
}

impl Mate {
    fn label(self) -> &'static str {
        match self {
            Self::R1 => "R1",
            Self::R2 => "R2",
        }
    }

    fn thread_name(self) -> &'static str {
        match self {
            Self::R1 => "pairasm-fastq-r1-reader",
            Self::R2 => "pairasm-fastq-r2-reader",
        }
    }
}

struct MateBatch {
    records: Vec<FastqRecord>,
    active_len: usize,
}

impl MateBatch {
    fn new() -> Self {
        Self {
            records: (0..IO_BATCH_SIZE).map(|_| FastqRecord::default()).collect(),
            active_len: 0,
        }
    }

    fn active_records(&self) -> &[FastqRecord] {
        &self.records[..self.active_len]
    }

    fn len(&self) -> usize {
        self.active_len
    }
}

enum MateBatchMessage {
    Batch(MateBatch),
    EndOfInput,
}

struct MateReaderWorker {
    mate: Mate,
    empty_batches: Option<SyncSender<MateBatch>>,
    filled_batches: Receiver<Result<MateBatchMessage>>,
    handle: Option<JoinHandle<()>>,
}

impl MateReaderWorker {
    fn spawn(mate: Mate, path: &str) -> Result<Self> {
        let (empty_tx, empty_rx) = mpsc::sync_channel::<MateBatch>(1);
        let (filled_tx, filled_rx) = mpsc::sync_channel::<Result<MateBatchMessage>>(1);
        let thread_path = path.to_owned();

        let handle = thread::Builder::new()
            .name(mate.thread_name().to_owned())
            .spawn(move || {
                let mut reader = match open_fastq_reader(&thread_path) {
                    Ok(reader) => reader,
                    Err(error) => {
                        let _ = filled_tx.send(Err(error));
                        return;
                    },
                };

                while let Ok(mut batch) = empty_rx.recv() {
                    match fill_mate_batch(&mut reader, &thread_path, &mut batch) {
                        Ok(PairReadStatus::Read) => {
                            if filled_tx.send(Ok(MateBatchMessage::Batch(batch))).is_err() {
                                return;
                            }
                        },
                        Ok(PairReadStatus::EndOfInput) => {
                            let _ = filled_tx.send(Ok(MateBatchMessage::EndOfInput));
                            return;
                        },
                        Err(error) => {
                            let _ = filled_tx.send(Err(error));
                            return;
                        },
                    }
                }
            })
            .wrap_err_with(|| {
                format!(
                    "failed to spawn {} FASTQ reader thread for {path}",
                    mate.label()
                )
            })?;

        let worker = Self {
            mate,
            empty_batches: Some(empty_tx),
            filled_batches: filled_rx,
            handle: Some(handle),
        };
        worker.send_empty_batch(MateBatch::new())?;
        Ok(worker)
    }

    fn send_empty_batch(&self, batch: MateBatch) -> Result<()> {
        let empty_batches = self.empty_batches.as_ref().ok_or_else(|| {
            eyre!(
                "{} FASTQ reader thread was cancelled before receiving a reusable batch",
                self.mate.label()
            )
        })?;

        empty_batches.send(batch).wrap_err_with(|| {
            format!(
                "{} FASTQ reader thread stopped before receiving a reusable batch",
                self.mate.label()
            )
        })
    }

    fn recv_filled_batch(&self) -> Result<MateBatchMessage> {
        self.filled_batches.recv().wrap_err_with(|| {
            format!(
                "{} FASTQ reader thread stopped before sending a batch",
                self.mate.label()
            )
        })?
    }

    fn cancel(&mut self) {
        self.empty_batches.take();
    }

    fn join(mut self) -> Result<()> {
        self.cancel();
        self.handle
            .take()
            .expect("FASTQ reader worker should be joined at most once")
            .join()
            .map_err(|_| eyre!("{} FASTQ reader thread panicked", self.mate.label()))
    }
}

struct PairedBatch {
    r1: MateBatch,
    r2: MateBatch,
}

impl PairedBatch {
    fn new(
        r1: MateBatch,
        r2: MateBatch,
        input_label: &str,
        complete_pairs_seen: u64,
    ) -> Result<Self> {
        if r1.len() != r2.len() {
            bail!(
                "paired FASTQ inputs have different record counts\nsource: {input_label}\ncomplete_pairs_seen: {complete_pairs_seen}\nhelp: pairasm expects R1 and R2 FASTQs to be in the same order and have the same number of records",
            );
        }

        Ok(Self { r1, r2 })
    }

    fn from_messages(
        r1: MateBatchMessage,
        r2: MateBatchMessage,
        input_label: &str,
        complete_pairs_seen: u64,
    ) -> Result<Option<Self>> {
        match (r1, r2) {
            (MateBatchMessage::Batch(r1), MateBatchMessage::Batch(r2)) => {
                Self::new(r1, r2, input_label, complete_pairs_seen).map(Some)
            },
            (MateBatchMessage::EndOfInput, MateBatchMessage::EndOfInput) => Ok(None),
            (MateBatchMessage::Batch(_), MateBatchMessage::EndOfInput)
            | (MateBatchMessage::EndOfInput, MateBatchMessage::Batch(_)) => {
                bail!(
                    "paired FASTQ inputs have different record counts\nsource: {input_label}\ncomplete_pairs_seen: {complete_pairs_seen}\nhelp: pairasm expects R1 and R2 FASTQs to be in the same order and have the same number of records",
                );
            },
        }
    }

    fn len(&self) -> usize {
        self.r1.len()
    }

    fn records(&self) -> impl Iterator<Item = (&FastqRecord, &FastqRecord)> + '_ {
        self.r1
            .active_records()
            .iter()
            .zip(self.r2.active_records())
    }

    fn r1_records(&self) -> &[FastqRecord] {
        self.r1.active_records()
    }

    fn r2_records(&self) -> &[FastqRecord] {
        self.r2.active_records()
    }

    fn into_inner(self) -> (MateBatch, MateBatch) {
        (self.r1, self.r2)
    }
}

struct PairedFastqReaders {
    r1: MateReaderWorker,
    r2: MateReaderWorker,
    input_label: String,
    complete_pairs_seen: u64,
}

impl PairedFastqReaders {
    fn spawn(input1: &str, input2: &str, input_label: String) -> Result<Self> {
        Ok(Self {
            r1: MateReaderWorker::spawn(Mate::R1, input1)?,
            r2: MateReaderWorker::spawn(Mate::R2, input2)?,
            input_label,
            complete_pairs_seen: 0,
        })
    }

    fn next_batch(&mut self) -> Result<Option<PairedBatch>> {
        let batch = PairedBatch::from_messages(
            self.r1.recv_filled_batch()?,
            self.r2.recv_filled_batch()?,
            &self.input_label,
            self.complete_pairs_seen,
        )?;

        if let Some(batch) = &batch {
            self.complete_pairs_seen += batch.len() as u64;
        }

        Ok(batch)
    }

    fn recycle(&mut self, batch: PairedBatch) -> Result<()> {
        let (mut r1, mut r2) = batch.into_inner();
        r1.active_len = 0;
        r2.active_len = 0;

        if let Err(error) = self.r1.send_empty_batch(r1) {
            self.cancel();
            return Err(error);
        }
        if let Err(error) = self.r2.send_empty_batch(r2) {
            self.cancel();
            return Err(error);
        }

        Ok(())
    }

    fn cancel(&mut self) {
        self.r1.cancel();
        self.r2.cancel();
    }

    fn finish(self) -> Result<()> {
        let r1_result = self.r1.join();
        let r2_result = self.r2.join();
        r1_result?;
        r2_result
    }
}

fn fill_mate_batch(
    reader: &mut FastqReader<Box<dyn BufRead>>,
    path: &str,
    batch: &mut MateBatch,
) -> Result<PairReadStatus> {
    batch.active_len = 0;

    while batch.active_len < batch.records.len() {
        let bytes_read = reader
            .read_record(&mut batch.records[batch.active_len])
            .wrap_err_with(|| format!("failed to read {path}"))?;

        if bytes_read == 0 {
            break;
        }

        batch.active_len += 1;
    }

    if batch.active_len == 0 {
        Ok(PairReadStatus::EndOfInput)
    } else {
        Ok(PairReadStatus::Read)
    }
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

    fn write_merged(&mut self, read: &OwnedSequenceRead) -> Result<()> {
        self.merged
            .write_record(&merged_record(read))
            .wrap_err_with(|| {
                format!(
                    "failed to write merged FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                    read.id()
                )
            })
    }

    fn write_unmerged_records(&mut self, r1: &FastqRecord, r2: &FastqRecord) -> Result<bool> {
        let Some(output) = &mut self.unmerged else {
            return Ok(false);
        };

        output.write_record(r1).wrap_err_with(|| {
            format!(
                "failed to write unmerged R1 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(r1.name().as_ref())
            )
        })?;
        output.write_record(r2).wrap_err_with(|| {
            format!(
                "failed to write unmerged R2 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(r2.name().as_ref())
            )
        })?;

        Ok(true)
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

#[derive(Clone, Copy)]
struct FastqMateKey<'key>(&'key [u8]);

impl<'key> FastqMateKey<'key> {
    fn from_record(record: &'key FastqRecord) -> Self {
        Self::from_header(record.name().as_ref())
    }

    fn from_header(header: &'key [u8]) -> Self {
        let first_token = header
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(header);

        let key = match first_token {
            [prefix @ .., b'/', b'1' | b'2'] => prefix,
            _ => first_token,
        };

        Self(key)
    }

    fn as_bytes(self) -> &'key [u8] {
        self.0
    }

    fn as_str(self) -> Result<&'key str> {
        str::from_utf8(self.0).wrap_err("FASTQ read name is not valid UTF-8")
    }
}

fn mate_keys_match(r1: &FastqRecord, r2: &FastqRecord) -> bool {
    FastqMateKey::from_record(r1).as_bytes() == FastqMateKey::from_record(r2).as_bytes()
}

struct FastqReadView<'read> {
    id: &'read str,
    seq: &'read str,
    qual: &'read str,
}

impl<'read> FastqReadView<'read> {
    fn from_record(record: &'read FastqRecord, id: &'read str) -> Result<Self> {
        Ok(Self {
            id,
            seq: str::from_utf8(record.sequence()).wrap_err("FASTQ sequence is not valid UTF-8")?,
            qual: str::from_utf8(record.quality_scores())
                .wrap_err("FASTQ quality string is not valid UTF-8")?,
        })
    }
}

impl SeqRecordView for FastqReadView<'_> {
    fn id(&self) -> &str {
        self.id
    }

    fn seq(&self) -> &str {
        self.seq
    }

    fn qual(&self) -> &str {
        self.qual
    }
}

fn merged_record(read: &OwnedSequenceRead) -> FastqRecord {
    FastqRecord::new(
        Definition::new(read.id(), ""),
        read.sequence_bytes(),
        read.quality_bytes(),
    )
}

struct MergeOrchestrator<'request> {
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
    fn new(request: &'request RunRequest) -> Result<Self> {
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

    fn run(mut self) -> Result<()> {
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

            if let Err(error) = self.inputs.recycle(batch) {
                self.inputs.cancel();
                return Err(error);
            }
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

enum PairReadStatus {
    Read,
    EndOfInput,
}

/// Run pair merging over two FASTQ inputs.
///
/// # Errors
///
/// Returns an error when input files cannot be read, paired inputs violate the ordering
/// contract too often, output cannot be written, or a non-biological assembly invariant fails.
pub fn run(request: &RunRequest) -> Result<()> {
    MergeOrchestrator::new(request)?.run()
}

#[cfg(test)]
mod tests {
    use super::{Definition, FastqMateKey, FastqRecord, MateBatch, MateBatchMessage, PairedBatch};

    #[test]
    fn mate_key_strips_slash_mate_suffix() {
        assert_eq!(
            FastqMateKey::from_header(b"read123/1").as_bytes(),
            b"read123"
        );
        assert_eq!(
            FastqMateKey::from_header(b"read123/2").as_bytes(),
            b"read123"
        );
        assert_eq!(FastqMateKey::from_header(b"read123").as_bytes(), b"read123");
    }

    #[test]
    fn mate_key_uses_first_whitespace_token() {
        assert_eq!(
            FastqMateKey::from_header(b"read123/1 instrument stuff").as_bytes(),
            b"read123"
        );
        assert_eq!(
            FastqMateKey::from_header(b"read123 comment").as_bytes(),
            b"read123"
        );
    }

    #[test]
    fn mate_key_can_be_extracted_from_fastq_record() {
        let record = FastqRecord::new(
            Definition::new("read123/1", "instrument stuff"),
            "AAAA",
            "IIII",
        );

        assert_eq!(FastqMateKey::from_record(&record).as_bytes(), b"read123");
    }

    #[test]
    fn paired_batch_from_two_eof_messages_finishes_input() {
        assert!(
            PairedBatch::from_messages(
                MateBatchMessage::EndOfInput,
                MateBatchMessage::EndOfInput,
                "test-input",
                42,
            )
            .expect("EOF messages should parse successfully")
            .is_none()
        );
    }

    #[test]
    fn paired_batch_rejects_different_active_lengths() {
        let Err(error) = PairedBatch::from_messages(
            MateBatchMessage::Batch(mate_batch_with_len(2)),
            MateBatchMessage::Batch(mate_batch_with_len(1)),
            "test-input",
            42,
        ) else {
            panic!("paired batch with different active lengths should fail");
        };

        assert!(
            error
                .to_string()
                .contains("paired FASTQ inputs have different record counts")
        );
    }

    #[test]
    fn paired_batch_accepts_equal_active_lengths() {
        let batch = PairedBatch::from_messages(
            MateBatchMessage::Batch(mate_batch_with_len(2)),
            MateBatchMessage::Batch(mate_batch_with_len(2)),
            "test-input",
            42,
        )
        .expect("equal active lengths should parse successfully")
        .expect("equal active lengths should produce a paired batch");

        assert_eq!(batch.len(), 2);
        assert_eq!(batch.records().count(), 2);
    }

    fn mate_batch_with_len(active_len: usize) -> MateBatch {
        let mut batch = MateBatch::new();
        batch.active_len = active_len;
        batch
    }
}
