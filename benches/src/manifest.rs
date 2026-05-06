use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use color_eyre::eyre::{Result, WrapErr, bail};
use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};

pub fn write_merged_manifest(merged_fastq: &Path, manifest_path: &Path) -> Result<usize> {
    let mut writer = BufWriter::new(File::create(manifest_path)?);
    writeln!(
        writer,
        "read_id\toutput_header\tmerged_len\tavg_qual\tmin_qual\tmax_qual\tsequence_hash\tquality_hash"
    )?;

    if !merged_fastq.exists() {
        writer.flush()?;
        return Ok(0);
    }

    let mut reader = open_fastq_reader(merged_fastq)?;
    let mut rows = 0usize;
    while let Some(record) = read_record(&mut reader, merged_fastq)? {
        let quality = quality_summary(&record.quality)?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            clean_tsv_field(&normalize_read_id(&record.header)),
            clean_tsv_field(&record.header),
            record.sequence.len(),
            quality.avg,
            quality.min,
            quality.max,
            stable_hash(record.sequence.as_bytes()),
            stable_hash(record.quality.as_bytes())
        )?;
        rows += 1;
    }

    writer.flush()?;
    Ok(rows)
}

fn open_fastq_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).wrap_err_with(|| format!("failed to open {}", path.display()))?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn read_record(reader: &mut dyn BufRead, path: &Path) -> Result<Option<FastqRecord>> {
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(None);
    }

    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    if reader.read_line(&mut sequence)? == 0
        || reader.read_line(&mut plus)? == 0
        || reader.read_line(&mut quality)? == 0
    {
        bail!("truncated FASTQ record in {}", path.display());
    }

    trim_line_end(&mut header);
    trim_line_end(&mut sequence);
    trim_line_end(&mut plus);
    trim_line_end(&mut quality);

    if !header.starts_with('@') {
        bail!(
            "FASTQ record header does not start with @ in {}",
            path.display()
        );
    }
    if !plus.starts_with('+') {
        bail!(
            "FASTQ record separator does not start with + in {}",
            path.display()
        );
    }
    if sequence.len() != quality.len() {
        bail!(
            "FASTQ sequence and quality lengths differ in {} for {}: seq_len={}, qual_len={}",
            path.display(),
            header,
            sequence.len(),
            quality.len()
        );
    }

    Ok(Some(FastqRecord {
        header,
        sequence,
        quality,
    }))
}

fn trim_line_end(line: &mut String) {
    while line.ends_with('\n') || line.ends_with('\r') {
        line.pop();
    }
}

fn normalize_read_id(header: &str) -> String {
    let token = header
        .strip_prefix('@')
        .unwrap_or(header)
        .split_whitespace()
        .next()
        .unwrap_or_default();

    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
        .to_owned()
}

fn clean_tsv_field(field: &str) -> String {
    field
        .chars()
        .map(|ch| {
            if matches!(ch, '\t' | '\r' | '\n') {
                ' '
            } else {
                ch
            }
        })
        .collect()
}

fn stable_hash(bytes: &[u8]) -> String {
    format!("sha256:{}", hex::encode(Sha256::digest(bytes)))
}

fn quality_summary(quality_ascii: &str) -> Result<QualitySummary> {
    let mut min = u8::MAX;
    let mut max = 0u8;
    let mut sum = 0u64;

    for quality in quality_ascii.bytes() {
        if quality < 33 {
            bail!("FASTQ quality byte {quality} is below Phred+33 ASCII range");
        }
        let score = quality - 33;
        min = min.min(score);
        max = max.max(score);
        sum += u64::from(score);
    }

    if quality_ascii.is_empty() {
        return Ok(QualitySummary {
            avg: 0.0,
            min: 0,
            max: 0,
        });
    }

    Ok(QualitySummary {
        avg: average_quality(sum, quality_ascii.len()),
        min,
        max,
    })
}

#[expect(
    clippy::cast_precision_loss,
    reason = "manifest quality means are descriptive benchmark summaries, not exact identifiers"
)]
fn average_quality(sum: u64, len: usize) -> f64 {
    sum as f64 / len as f64
}

struct FastqRecord {
    header: String,
    sequence: String,
    quality: String,
}

struct QualitySummary {
    avg: f64,
    min: u8,
    max: u8,
}

#[cfg(test)]
mod tests {
    use super::{normalize_read_id, quality_summary, stable_hash};

    #[test]
    fn normalizes_pair_id_like_pairasm_cli() {
        assert_eq!(normalize_read_id("@read/1 extra metadata"), "read");
        assert_eq!(normalize_read_id("@read/2"), "read");
        assert_eq!(normalize_read_id("@read stuff/1"), "read");
    }

    #[test]
    fn quality_summary_reports_phred_scores() {
        let summary = quality_summary("!+5").expect("quality summary should parse");
        assert_eq!(summary.min, 0);
        assert_eq!(summary.max, 20);
        assert!((summary.avg - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn stable_hash_is_prefixed_sha256() {
        assert_eq!(
            stable_hash(b"ACGT"),
            "sha256:1dff3e84fe7877e0673b69bbddcf40124e396e3f9943dd890c91b6a09adb9af0"
        );
    }
}
