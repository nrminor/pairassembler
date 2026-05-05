use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use color_eyre::eyre::Result;
use flate2::{Compression, read::GzDecoder, write::GzEncoder};

pub fn validate_gzip(path: &Path) -> Result<()> {
    let mut decoder = GzDecoder::new(File::open(path)?);
    io::copy(&mut decoder, &mut io::sink())?;
    Ok(())
}

pub fn fastq_record_count(path: &Path) -> Result<usize> {
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|extension| extension == "gz") {
        Box::new(BufReader::new(GzDecoder::new(File::open(path)?)))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };

    let mut lines = 0usize;
    for line in reader.lines() {
        let _ = line?;
        lines += 1;
    }
    Ok(lines / 4)
}

pub fn file_size(path: &Path) -> Result<u64> {
    Ok(path.metadata()?.len())
}

pub fn write_first_fastq_records(src: &Path, dst: &Path, read_pairs: usize) -> Result<()> {
    let input = File::open(src)?;
    let decoder = GzDecoder::new(input);
    let mut reader = BufReader::new(decoder);
    let output = File::create(dst)?;
    let encoder = GzEncoder::new(output, Compression::default());
    let mut writer = BufWriter::new(encoder);
    let mut line = String::new();
    let mut lines_written = 0usize;
    let max_lines = read_pairs.saturating_mul(4);

    while lines_written < max_lines {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        writer.write_all(line.as_bytes())?;
        lines_written += 1;
    }
    writer.flush()?;
    Ok(())
}
