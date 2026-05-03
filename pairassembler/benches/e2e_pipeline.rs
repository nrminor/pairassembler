use std::{
    error::Error,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
    str,
    time::{Duration, Instant},
};

use criterion::{Criterion, criterion_group, criterion_main};
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use libpairassembly::{
    Assembler, OverlapParams, OverlapSearch, OverlapValidator, OwnedSequenceRead, PairInput,
    SequenceRead,
};
use noodles::fastq::{Record as FastqRecord, io::Writer as FastqWriter, record::Definition};
use pairassembler::{RunRequest, RunSettings, cli::UiPolicy, merging, progress::ProgressMode};
use tempfile::TempDir;

const DEFAULT_E2E_PAIRS: usize = 10_000;

struct FastqPair {
    id: String,
    r1_sequence: String,
    r2_sequence: String,
    quality: String,
}

fn e2e_plain_pipeline(c: &mut Criterion) {
    let temp = tempdir_or_panic();
    let pair_count = e2e_pair_count();
    let (r1, r2) = write_fastq_pair_files(temp.path(), "plain", pair_count)
        .unwrap_or_else(|error| panic!("failed to write benchmark FASTQ inputs: {error}"));
    let merged = temp.path().join("merged.fastq");
    let label = format!("e2e_plain_mergeable_{pair_count}");
    print_phase_report(&label, &r1, &r2, &merged);

    c.bench_function(&label, |b| {
        b.iter(|| run_pipeline(&r1, &r2, &merged));
    });
}

fn e2e_gzip_pipeline(c: &mut Criterion) {
    let temp = tempdir_or_panic();
    let pair_count = e2e_pair_count();
    let (r1, r2) = write_gzip_fastq_pair_files(temp.path(), "gzip", pair_count)
        .unwrap_or_else(|error| panic!("failed to write gzip benchmark inputs: {error}"));
    let merged = temp.path().join("merged.fastq.gz");
    let label = format!("e2e_gzip_mergeable_{pair_count}");
    print_phase_report(&label, &r1, &r2, &merged);

    c.bench_function(&label, |b| {
        b.iter(|| run_pipeline(&r1, &r2, &merged));
    });
}

fn e2e_pair_count() -> usize {
    std::env::var("PAIRASM_E2E_PAIRS")
        .ok()
        .and_then(|raw| raw.parse::<usize>().ok())
        .filter(|count| *count > 0)
        .unwrap_or(DEFAULT_E2E_PAIRS)
}

fn run_pipeline(r1: &Path, r2: &Path, merged: &Path) {
    let request = RunRequest {
        input1: path_to_string(r1),
        input2: path_to_string(r2),
        output_file: Some(path_to_string(merged)),
        unmerged_output: None,
        summary: None,
        progress_every: 0,
        ui: UiPolicy {
            log_level: None,
            show_summary: false,
            progress_mode: ProgressMode::Off,
        },
        settings: RunSettings::new(
            OverlapParams::default(),
            OverlapValidator::default(),
            false,
            3,
        ),
    };

    if let Err(error) = merging::run(&request) {
        panic!("unexpected e2e benchmark pipeline failure: {error}");
    }
}

type BenchResult<T> = Result<T, Box<dyn Error>>;

enum BenchFastqOutput {
    Plain(FastqWriter<Box<dyn Write>>),
    Gzip(FastqWriter<GzEncoder<Box<dyn Write>>>),
}

impl BenchFastqOutput {
    fn new(path: &Path) -> BenchResult<Self> {
        let file = File::create(path)?;
        let writer: Box<dyn Write> = Box::new(BufWriter::new(file));
        if is_gzip_path(path) {
            Ok(Self::Gzip(FastqWriter::new(GzEncoder::new(
                writer,
                Compression::default(),
            ))))
        } else {
            Ok(Self::Plain(FastqWriter::new(writer)))
        }
    }

    fn write_record(&mut self, record: &FastqRecord) -> BenchResult<()> {
        match self {
            Self::Plain(writer) => writer.write_record(record)?,
            Self::Gzip(writer) => writer.write_record(record)?,
        }
        Ok(())
    }

    fn finish(self) -> BenchResult<()> {
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

#[derive(Clone, Copy, Debug, Default)]
struct PhaseTimings {
    setup: Duration,
    read_decode_parse: Duration,
    assemble: Duration,
    encode_write: Duration,
    finish_flush: Duration,
    pairs_seen: usize,
    pairs_merged: usize,
}

impl PhaseTimings {
    fn total(self) -> Duration {
        self.setup + self.read_decode_parse + self.assemble + self.encode_write + self.finish_flush
    }
}

fn print_phase_report(label: &str, r1: &Path, r2: &Path, merged: &Path) {
    let timings = measure_pipeline_phases(r1, r2, merged)
        .unwrap_or_else(|error| panic!("failed to measure benchmark phases for {label}: {error}"));
    let total = timings.total();
    eprintln!("\nphase report: {label}");
    eprintln!("  pairs seen:     {}", timings.pairs_seen);
    eprintln!("  pairs merged:   {}", timings.pairs_merged);
    eprintln!("  setup:          {}", phase_line(timings.setup, total));
    eprintln!(
        "  read+decode+parse: {}",
        phase_line(timings.read_decode_parse, total)
    );
    eprintln!("  assemble:       {}", phase_line(timings.assemble, total));
    eprintln!(
        "  encode+write:   {}",
        phase_line(timings.encode_write, total)
    );
    eprintln!(
        "  finish+flush:   {}",
        phase_line(timings.finish_flush, total)
    );
    eprintln!(
        "  total measured: {:.3} ms\n",
        total.as_secs_f64() * 1_000.0
    );
}

fn phase_line(phase: Duration, total: Duration) -> String {
    let phase_ms = phase.as_secs_f64() * 1_000.0;
    let percent = if total.is_zero() {
        0.0
    } else {
        phase.as_secs_f64() / total.as_secs_f64() * 100.0
    };
    format!("{phase_ms:>8.3} ms ({percent:>5.1}%)")
}

fn measure_pipeline_phases(r1: &Path, r2: &Path, merged: &Path) -> BenchResult<PhaseTimings> {
    let mut timings = PhaseTimings::default();

    let phase_start = Instant::now();
    let mut fastq_reader1 = open_fastq_reader(r1)?;
    let mut fastq_reader2 = open_fastq_reader(r2)?;
    let mut records1 = fastq_reader1.records();
    let mut records2 = fastq_reader2.records();
    let mut output = BenchFastqOutput::new(merged)?;
    let assembler = Assembler::builder().build()?;
    timings.setup += phase_start.elapsed();

    loop {
        let phase_start = Instant::now();
        let next1 = records1.next();
        let next2 = records2.next();
        timings.read_decode_parse += phase_start.elapsed();

        let Some((fwd, rev)) = next_pair(next1, next2)? else {
            break;
        };
        timings.pairs_seen = timings.pairs_seen.saturating_add(1);

        let phase_start = Instant::now();
        let merged_read = assemble_record_pair(&assembler, &fwd, &rev)?;
        timings.assemble += phase_start.elapsed();

        if let Some(read) = merged_read {
            let record = merged_record(&read);
            let phase_start = Instant::now();
            output.write_record(&record)?;
            timings.encode_write += phase_start.elapsed();
            timings.pairs_merged = timings.pairs_merged.saturating_add(1);
        }
    }

    let phase_start = Instant::now();
    output.finish()?;
    timings.finish_flush += phase_start.elapsed();

    Ok(timings)
}

fn open_fastq_reader(path: &Path) -> BenchResult<noodles::fastq::Reader<Box<dyn BufRead>>> {
    let read_buffer = BufReader::new(File::open(path)?);
    let reader: Box<dyn BufRead> = if is_gzip_path(path) {
        Box::new(BufReader::new(GzDecoder::new(read_buffer)))
    } else {
        Box::new(read_buffer)
    };
    Ok(noodles::fastq::Reader::new(reader))
}

fn next_pair(
    next1: Option<io::Result<FastqRecord>>,
    next2: Option<io::Result<FastqRecord>>,
) -> BenchResult<Option<(FastqRecord, FastqRecord)>> {
    match (next1, next2) {
        (Some(Ok(fwd)), Some(Ok(rev))) => Ok(Some((fwd, rev))),
        (None, None) => Ok(None),
        (Some(Err(error)), _) | (_, Some(Err(error))) => Err(Box::new(error)),
        (Some(Ok(_)), None) | (None, Some(Ok(_))) => {
            Err("benchmark FASTQ inputs have different record counts".into())
        },
    }
}

fn assemble_record_pair(
    assembler: &Assembler,
    fwd: &FastqRecord,
    rev: &FastqRecord,
) -> BenchResult<Option<OwnedSequenceRead>> {
    let pair_key = str::from_utf8(mate_key(fwd.name().as_ref()))?;
    let read1 = sequence_read_from_record(fwd, pair_key)?;
    let read2 = sequence_read_from_record(rev, pair_key)?;
    let pair_input = PairInput::new(read1, read2);

    let overlap = match assembler.on_pair(&pair_input)?.find_overlap()? {
        OverlapSearch::Found(ctx) => ctx,
        OverlapSearch::NoOverlap(_) => return Ok(None),
    };
    let Ok(validated) = overlap.validate() else {
        return Ok(None);
    };
    Ok(Some(validated.merge()?.correct()?.into_owned_read()?))
}

fn sequence_read_from_record<'a>(
    record: &'a FastqRecord,
    pair_id: &'a str,
) -> BenchResult<SequenceRead<'a>> {
    let seq = str::from_utf8(record.sequence())?;
    let qual = str::from_utf8(record.quality_scores())?;
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

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|extension| extension == "gz")
}

fn path_to_string(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn tempdir_or_panic() -> TempDir {
    TempDir::new().unwrap_or_else(|error| panic!("failed to create benchmark tempdir: {error}"))
}

fn write_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pair_count: usize,
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq"));
    let r2 = directory.join(format!("{stem}_R2.fastq"));
    let mut r1_writer = File::create(&r1)?;
    let mut r2_writer = File::create(&r2)?;
    write_pairs(&mut r1_writer, &mut r2_writer, pair_count)?;
    Ok((r1, r2))
}

fn write_gzip_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pair_count: usize,
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq.gz"));
    let r2 = directory.join(format!("{stem}_R2.fastq.gz"));
    let mut r1_writer = GzEncoder::new(File::create(&r1)?, Compression::default());
    let mut r2_writer = GzEncoder::new(File::create(&r2)?, Compression::default());
    write_pairs(&mut r1_writer, &mut r2_writer, pair_count)?;
    let _r1_file = r1_writer.finish()?;
    let _r2_file = r2_writer.finish()?;
    Ok((r1, r2))
}

fn write_pairs(
    r1_writer: &mut impl Write,
    r2_writer: &mut impl Write,
    pair_count: usize,
) -> io::Result<()> {
    for index in 0..pair_count {
        let pair = mergeable_pair(index);
        write_record(
            r1_writer,
            &format!("{}/1", pair.id),
            &pair.r1_sequence,
            &pair.quality,
        )?;
        write_record(
            r2_writer,
            &format!("{}/2", pair.id),
            &pair.r2_sequence,
            &pair.quality,
        )?;
    }
    Ok(())
}

fn mergeable_pair(index: usize) -> FastqPair {
    FastqPair::new(
        &format!("mergeable-{index:08}"),
        "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
        "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
    )
}

fn write_record(
    writer: &mut impl Write,
    id: &str,
    sequence: &str,
    quality: &str,
) -> io::Result<()> {
    writeln!(writer, "@{id}\n{sequence}\n+\n{quality}")
}

impl FastqPair {
    fn new(id: &str, r1_sequence: &str, r2_sequence: &str) -> Self {
        Self {
            id: id.to_owned(),
            r1_sequence: r1_sequence.to_owned(),
            r2_sequence: r2_sequence.to_owned(),
            quality: "I".repeat(r1_sequence.len()),
        }
    }
}

criterion_group!(benches, e2e_plain_pipeline, e2e_gzip_pipeline);
criterion_main!(benches);
