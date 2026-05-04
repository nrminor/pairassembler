use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    mem,
    path::Path,
    str,
    time::Instant,
};

use color_eyre::eyre::{Result, WrapErr, bail};
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

    fn write_unmerged_pair(&mut self, pair: &PairInput<FastqRecord>) -> Result<bool> {
        let Some(output) = &mut self.unmerged else {
            return Ok(false);
        };

        output.write_record(&pair.r1).wrap_err_with(|| {
            format!(
                "failed to write unmerged R1 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(pair.r1.name().as_ref())
            )
        })?;
        output.write_record(&pair.r2).wrap_err_with(|| {
            format!(
                "failed to write unmerged R2 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(pair.r2.name().as_ref())
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
    reader1: FastqReader<Box<dyn BufRead>>,
    reader2: FastqReader<Box<dyn BufRead>>,
    outputs: OutputHandles,
    assembler: Assembler,
    no_correct: bool,
    context: RunContext,
    stats: AssemblyStats,
    progress: ProgressReporter,
    started_at: Instant,
    input_batch: Vec<PairInput<FastqRecord>>,
    outcomes: Vec<Result<PairAssemblyOutcome>>,
}

impl<'request> MergeOrchestrator<'request> {
    fn new(request: &'request RunRequest) -> Result<Self> {
        let input_batch = (0..IO_BATCH_SIZE)
            .map(|_| PairInput::new(FastqRecord::default(), FastqRecord::default()))
            .collect();

        Ok(Self {
            request,
            reader1: open_fastq_reader(&request.input1)?,
            reader2: open_fastq_reader(&request.input2)?,
            outputs: OutputHandles::new(request)?,
            assembler: Assembler::builder()
                .with_overlap_params(request.settings.overlap_settings)
                .with_validator(request.settings.validation_settings)
                .build()?,
            no_correct: request.settings.no_correct,
            context: RunContext::from_request(request),
            stats: AssemblyStats::new(!request.settings.no_correct),
            progress: ProgressReporter::new(request.ui.progress_mode, request.progress_every),
            started_at: Instant::now(),
            input_batch,
            outcomes: Vec::with_capacity(IO_BATCH_SIZE),
        })
    }

    fn run(mut self) -> Result<()> {
        loop {
            let active_len = self.fill_input_batch()?;

            if active_len == 0 {
                break;
            }

            self.assemble_active_batch(active_len);
            self.write_outcomes(active_len)?;
        }

        self.finish()
    }

    fn fill_input_batch(&mut self) -> Result<usize> {
        let mut active_len = 0;

        while active_len < self.input_batch.len() {
            match self.read_pair_into_slot(active_len)? {
                PairReadStatus::Read => {},
                PairReadStatus::EndOfInput => break,
            }

            self.record_pair_seen(active_len);

            if self.record_mate_mismatch_if_needed(active_len)? {
                self.progress.maybe_report(&self.stats);
                continue;
            }

            active_len += 1;
        }

        Ok(active_len)
    }

    fn assemble_active_batch(&mut self, active_len: usize) {
        let assembler = &self.assembler;
        let no_correct = self.no_correct;
        let mut outcomes = mem::take(&mut self.outcomes);

        self.input_batch[..active_len]
            .par_iter()
            .with_min_len(MIN_PAIRS_PER_PARALLEL_TASK)
            .map(|pair| {
                let pair_id = FastqMateKey::from_record(&pair.r1).as_str()?;
                let pair_input = PairInput::new(
                    FastqReadView::from_record(&pair.r1, pair_id)?,
                    FastqReadView::from_record(&pair.r2, pair_id)?,
                );

                let overlap = match assembler.on_pair(&pair_input)?.find_overlap()? {
                    OverlapSearch::Found(overlap) => overlap,
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
                let read = if no_correct {
                    merged.into_owned_read()?
                } else {
                    merged.correct()?.into_owned_read()?
                };

                Ok(PairAssemblyOutcome::Merged(read))
            })
            .collect_into_vec(&mut outcomes);

        self.outcomes = outcomes;
    }

    fn read_pair_into_slot(&mut self, index: usize) -> Result<PairReadStatus> {
        let pair = &mut self.input_batch[index];
        let bytes_read1 = self
            .reader1
            .read_record(&mut pair.r1)
            .wrap_err_with(|| format!("failed to read {}", self.request.input1))?;
        let bytes_read2 = self
            .reader2
            .read_record(&mut pair.r2)
            .wrap_err_with(|| format!("failed to read {}", self.request.input2))?;

        match (bytes_read1, bytes_read2) {
            (0, 0) => Ok(PairReadStatus::EndOfInput),
            (0, _) | (_, 0) => {
                bail!(
                    "paired FASTQ inputs have different record counts\nsource: {}\ncomplete_pairs_seen: {}\nhelp: pairasm expects R1 and R2 FASTQs to be in the same order and have the same number of records",
                    self.context.input_label(),
                    self.stats.pairs_seen,
                );
            },
            (_, _) => Ok(PairReadStatus::Read),
        }
    }

    fn record_pair_seen(&mut self, index: usize) {
        let pair = &self.input_batch[index];
        self.stats
            .record_pair_seen(pair.r1.sequence().len(), pair.r2.sequence().len());
    }

    fn record_mate_mismatch_if_needed(&mut self, index: usize) -> Result<bool> {
        let pair = &self.input_batch[index];
        if FastqMateKey::from_record(&pair.r1).as_bytes()
            == FastqMateKey::from_record(&pair.r2).as_bytes()
        {
            return Ok(false);
        }

        self.stats.record_mate_id_mismatch();
        if self.stats.mate_id_mismatches > self.request.settings.max_mate_id_mismatches {
            bail!(
                "paired FASTQ inputs appear to be in different orders\nsource: {}\nmate_id_mismatches: {}\nmax_mate_id_mismatches: {}\nlast_r1_header: {}\nlast_r2_header: {}\npairs_seen_before_failure: {}\nhelp: pairasm expects R1 and R2 FASTQs to be sorted in the same pair order; repair or re-sort pairing before running pairasm",
                self.context.input_label(),
                self.stats.mate_id_mismatches,
                self.request.settings.max_mate_id_mismatches,
                String::from_utf8_lossy(pair.r1.name().as_ref()),
                String::from_utf8_lossy(pair.r2.name().as_ref()),
                self.stats.pairs_seen,
            );
        }

        Ok(true)
    }

    fn write_outcomes(&mut self, active_len: usize) -> Result<()> {
        for (pair, outcome) in self.input_batch[..active_len]
            .iter()
            .zip(self.outcomes.drain(..))
        {
            match outcome? {
                PairAssemblyOutcome::Merged(read) => {
                    self.outputs.write_merged(&read)?;
                    self.stats.record_merged(read.sequence_bytes().len());
                },
                PairAssemblyOutcome::Unmerged(reason) => {
                    let was_written = self.outputs.write_unmerged_pair(pair)?;
                    self.stats.record_unmerged(reason, was_written);
                },
            }

            self.progress.maybe_report(&self.stats);
        }

        Ok(())
    }

    fn finish(self) -> Result<()> {
        self.progress.finish();
        self.outputs.finish()?;

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
    use super::{Definition, FastqMateKey, FastqRecord};

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
}
