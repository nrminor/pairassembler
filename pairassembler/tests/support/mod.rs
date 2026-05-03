use std::{fs::File, io, io::Write, path::Path};

use flate2::{Compression, read::GzDecoder, write::GzEncoder};

#[derive(Clone, Copy, Debug)]
pub enum PairKind {
    Mergeable,
    NoOverlap,
    ValidationRejected,
}

#[derive(Clone, Debug)]
pub struct FastqPair {
    pub id: String,
    pub r1_sequence: String,
    pub r2_sequence: String,
    pub quality: String,
}

pub fn mixed_pairs() -> Vec<FastqPair> {
    vec![
        pair("mergeable-0001", PairKind::Mergeable),
        pair("no-overlap-0001", PairKind::NoOverlap),
        pair("validation-rejected-0001", PairKind::ValidationRejected),
    ]
}

pub fn many_pairs(count: usize, kind: PairKind) -> Vec<FastqPair> {
    (0..count)
        .map(|index| {
            let id = format!("{kind:?}-{index:08}");
            pair(&id, kind)
        })
        .collect()
}

pub fn write_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pairs: &[FastqPair],
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq"));
    let r2 = directory.join(format!("{stem}_R2.fastq"));
    write_fastq_pair_paths(&r1, &r2, pairs)?;
    Ok((r1, r2))
}

pub fn write_gzip_fastq_pair_files(
    directory: &Path,
    stem: &str,
    pairs: &[FastqPair],
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let r1 = directory.join(format!("{stem}_R1.fastq.gz"));
    let r2 = directory.join(format!("{stem}_R2.fastq.gz"));
    write_gzip_fastq_pair_paths(&r1, &r2, pairs)?;
    Ok((r1, r2))
}

pub fn write_fastq_pair_paths(r1: &Path, r2: &Path, pairs: &[FastqPair]) -> io::Result<()> {
    let mut r1_writer = File::create(r1)?;
    let mut r2_writer = File::create(r2)?;
    write_pairs(&mut r1_writer, &mut r2_writer, pairs)
}

pub fn write_gzip_fastq_pair_paths(r1: &Path, r2: &Path, pairs: &[FastqPair]) -> io::Result<()> {
    let r1_file = File::create(r1)?;
    let r2_file = File::create(r2)?;
    let mut r1_writer = GzEncoder::new(r1_file, Compression::default());
    let mut r2_writer = GzEncoder::new(r2_file, Compression::default());
    write_pairs(&mut r1_writer, &mut r2_writer, pairs)?;
    let _r1_file = r1_writer.finish()?;
    let _r2_file = r2_writer.finish()?;
    Ok(())
}

pub fn count_fastq_records(path: &Path) -> io::Result<usize> {
    let contents = if path.extension().is_some_and(|extension| extension == "gz") {
        let mut decoder = GzDecoder::new(File::open(path)?);
        let mut contents = String::new();
        std::io::Read::read_to_string(&mut decoder, &mut contents)?;
        contents
    } else {
        std::fs::read_to_string(path)?
    };

    Ok(contents.lines().count() / 4)
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
    writeln!(writer, "@{id}")?;
    writeln!(writer, "{sequence}")?;
    writeln!(writer, "+")?;
    writeln!(writer, "{quality}")?;
    Ok(())
}

fn pair(id: &str, kind: PairKind) -> FastqPair {
    match kind {
        PairKind::Mergeable => FastqPair::new(
            id,
            "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
            "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
        ),
        PairKind::NoOverlap => FastqPair::homopolymer(id, 'A', 'C', 56),
        PairKind::ValidationRejected => FastqPair::homopolymer(id, 'A', 'T', 56),
    }
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

    fn homopolymer(id: &str, r1_base: char, r2_base: char, len: usize) -> Self {
        Self {
            id: id.to_owned(),
            r1_sequence: r1_base.to_string().repeat(len),
            r2_sequence: r2_base.to_string().repeat(len),
            quality: "I".repeat(len),
        }
    }
}
