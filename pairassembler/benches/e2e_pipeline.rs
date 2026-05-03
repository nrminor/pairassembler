use std::{fs::File, io, io::Write, path::Path};

use criterion::{Criterion, criterion_group, criterion_main};
use flate2::{Compression, write::GzEncoder};
use libpairassembly::{OverlapParams, OverlapValidator};
use pairassembler::{RunRequest, RunSettings, cli::UiPolicy, merging, progress::ProgressMode};
use tempfile::TempDir;

struct FastqPair {
    id: String,
    r1_sequence: String,
    r2_sequence: String,
    quality: String,
}

fn e2e_plain_pipeline(c: &mut Criterion) {
    let temp = tempdir_or_panic();
    let pairs = many_mergeable_pairs(1_000);
    let (r1, r2) = write_fastq_pair_files(temp.path(), "plain", &pairs)
        .unwrap_or_else(|error| panic!("failed to write benchmark FASTQ inputs: {error}"));
    let merged = temp.path().join("merged.fastq");

    c.bench_function("e2e_plain_mergeable_1000", |b| {
        b.iter(|| run_pipeline(&r1, &r2, &merged));
    });
}

fn e2e_gzip_pipeline(c: &mut Criterion) {
    let temp = tempdir_or_panic();
    let pairs = many_mergeable_pairs(1_000);
    let (r1, r2) = write_gzip_fastq_pair_files(temp.path(), "gzip", &pairs)
        .unwrap_or_else(|error| panic!("failed to write gzip benchmark inputs: {error}"));
    let merged = temp.path().join("merged.fastq.gz");

    c.bench_function("e2e_gzip_mergeable_1000", |b| {
        b.iter(|| run_pipeline(&r1, &r2, &merged));
    });
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

fn path_to_string(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn tempdir_or_panic() -> TempDir {
    TempDir::new().unwrap_or_else(|error| panic!("failed to create benchmark tempdir: {error}"))
}

fn many_mergeable_pairs(count: usize) -> Vec<FastqPair> {
    (0..count)
        .map(|index| {
            FastqPair::new(
                &format!("mergeable-{index:08}"),
                "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
                "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
            )
        })
        .collect()
}

fn write_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pairs: &[FastqPair],
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq"));
    let r2 = directory.join(format!("{stem}_R2.fastq"));
    let mut r1_writer = File::create(&r1)?;
    let mut r2_writer = File::create(&r2)?;
    write_pairs(&mut r1_writer, &mut r2_writer, pairs)?;
    Ok((r1, r2))
}

fn write_gzip_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pairs: &[FastqPair],
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq.gz"));
    let r2 = directory.join(format!("{stem}_R2.fastq.gz"));
    let mut r1_writer = GzEncoder::new(File::create(&r1)?, Compression::default());
    let mut r2_writer = GzEncoder::new(File::create(&r2)?, Compression::default());
    write_pairs(&mut r1_writer, &mut r2_writer, pairs)?;
    let _r1_file = r1_writer.finish()?;
    let _r2_file = r2_writer.finish()?;
    Ok((r1, r2))
}

fn write_pairs(
    r1_writer: &mut impl Write,
    r2_writer: &mut impl Write,
    pairs: &[FastqPair],
) -> io::Result<()> {
    for pair in pairs {
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
