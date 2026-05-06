use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    sync::mpsc::{self, Receiver, SyncSender},
    thread::{self, JoinHandle},
};

use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use flate2::read::GzDecoder;
use noodles::fastq::{Reader as FastqReader, Record as FastqRecord};

use super::{IO_BATCH_SIZE, is_gzip_path};

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

pub(super) struct MateBatch {
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

pub(super) enum MateBatchMessage {
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

pub(super) struct PairedBatch {
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

    pub(super) fn len(&self) -> usize {
        self.r1.len()
    }

    pub(super) fn records(&self) -> impl Iterator<Item = (&FastqRecord, &FastqRecord)> + '_ {
        self.r1
            .active_records()
            .iter()
            .zip(self.r2.active_records())
    }

    pub(super) fn r1_records(&self) -> &[FastqRecord] {
        self.r1.active_records()
    }

    pub(super) fn r2_records(&self) -> &[FastqRecord] {
        self.r2.active_records()
    }

    fn into_inner(self) -> (MateBatch, MateBatch) {
        (self.r1, self.r2)
    }
}

pub(super) struct PairedFastqReaders {
    r1: MateReaderWorker,
    r2: MateReaderWorker,
    input_label: String,
    complete_pairs_seen: u64,
}

impl PairedFastqReaders {
    pub(super) fn spawn(input1: &str, input2: &str, input_label: String) -> Result<Self> {
        Ok(Self {
            r1: MateReaderWorker::spawn(Mate::R1, input1)?,
            r2: MateReaderWorker::spawn(Mate::R2, input2)?,
            input_label,
            complete_pairs_seen: 0,
        })
    }

    pub(super) fn next_batch(&mut self) -> Result<Option<PairedBatch>> {
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

    pub(super) fn recycle(&mut self, batch: PairedBatch) -> Result<()> {
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

    pub(super) fn cancel(&mut self) {
        self.r1.cancel();
        self.r2.cancel();
    }

    pub(super) fn finish(self) -> Result<()> {
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

enum PairReadStatus {
    Read,
    EndOfInput,
}

#[cfg(test)]
mod tests {
    use super::{MateBatch, MateBatchMessage, PairedBatch};

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
