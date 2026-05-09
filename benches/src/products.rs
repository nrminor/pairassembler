use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use color_eyre::eyre::{Result, WrapErr, bail};
use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};

pub struct MergedProduct {
    pub read_id: String,
    pub output_header: String,
    pub merged_len: usize,
    pub avg_qual: f64,
    pub min_qual: u8,
    pub max_qual: u8,
    pub sequence_hash: String,
    pub quality_hash: String,
}

pub fn read_merged_products(merged_fastq: &Path) -> Result<Vec<MergedProduct>> {
    if !merged_fastq.exists() {
        bail!("merged FASTQ does not exist: {}", merged_fastq.display());
    }

    let mut reader = open_fastq_reader(merged_fastq)?;
    let mut products = Vec::new();
    let mut read_ids = HashSet::new();
    let mut record = FastqRecord::default();
    while read_record(&mut reader, merged_fastq, &mut record)? {
        let quality = quality_summary(&record.quality)?;
        let read_id = normalize_read_id(&record.header);
        if !read_ids.insert(read_id.clone()) {
            bail!(
                "duplicate normalized merged read ID {read_id:?} in {}",
                merged_fastq.display()
            );
        }
        products.push(MergedProduct {
            read_id,
            output_header: record.header.clone(),
            merged_len: record.sequence.len(),
            avg_qual: quality.avg,
            min_qual: quality.min,
            max_qual: quality.max,
            sequence_hash: stable_hash(record.sequence.as_bytes()),
            quality_hash: stable_hash(record.quality.as_bytes()),
        });
    }

    Ok(products)
}

fn open_fastq_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).wrap_err_with(|| format!("failed to open {}", path.display()))?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn read_record(reader: &mut dyn BufRead, path: &Path, record: &mut FastqRecord) -> Result<bool> {
    record.header.clear();
    record.sequence.clear();
    record.plus.clear();
    record.quality.clear();

    if reader.read_line(&mut record.header)? == 0 {
        return Ok(false);
    }

    if reader.read_line(&mut record.sequence)? == 0
        || reader.read_line(&mut record.plus)? == 0
        || reader.read_line(&mut record.quality)? == 0
    {
        bail!("truncated FASTQ record in {}", path.display());
    }

    trim_line_end(&mut record.header);
    trim_line_end(&mut record.sequence);
    trim_line_end(&mut record.plus);
    trim_line_end(&mut record.quality);

    if !record.header.starts_with('@') {
        bail!(
            "FASTQ record header does not start with @ in {}",
            path.display()
        );
    }
    if !record.plus.starts_with('+') {
        bail!(
            "FASTQ record separator does not start with + in {}",
            path.display()
        );
    }
    if record.sequence.len() != record.quality.len() {
        bail!(
            "FASTQ sequence and quality lengths differ in {} for {}: seq_len={}, qual_len={}",
            path.display(),
            record.header,
            record.sequence.len(),
            record.quality.len()
        );
    }

    Ok(true)
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

fn stable_hash(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut encoded = [0u8; 64];
    hex::encode_to_slice(digest, &mut encoded).expect("SHA-256 digest hex encoding should fit");

    let mut hash = String::with_capacity("sha256:".len() + encoded.len());
    hash.push_str("sha256:");
    hash.push_str(std::str::from_utf8(&encoded).expect("hex output should be valid ASCII"));
    hash
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
    reason = "merged-product quality means are descriptive benchmark summaries, not exact identifiers"
)]
fn average_quality(sum: u64, len: usize) -> f64 {
    sum as f64 / len as f64
}

#[derive(Default)]
struct FastqRecord {
    header: String,
    sequence: String,
    plus: String,
    quality: String,
}

struct QualitySummary {
    avg: f64,
    min: u8,
    max: u8,
}

#[cfg(test)]
mod tests {
    use std::{fs, io::Write, path::PathBuf};

    use flate2::{Compression, write::GzEncoder};
    use uuid::Uuid;

    use super::{normalize_read_id, quality_summary, read_merged_products, stable_hash};

    #[test]
    fn normalizes_pair_id_like_pairasm_cli() {
        for (header, expected) in [
            ("@read/1 extra metadata", "read"),
            ("@read/2", "read"),
            ("@read stuff/1", "read"),
            ("@read/10 extra metadata", "read/10"),
            ("read/1 extra metadata", "read"),
            ("read/2", "read"),
            ("read/10", "read/10"),
            ("@read/1/2", "read/1"),
        ] {
            assert_eq!(normalize_read_id(header), expected, "header={header:?}");
        }
    }

    #[test]
    fn reads_plain_fastq_products() {
        let path = write_temp_fastq("@read/1 extra metadata\nACGT\n+\n!+5I\n@other/10\nA\n+\nI\n");

        let products = read_merged_products(&path).expect("FASTQ should parse");

        assert_eq!(products.len(), 2);
        assert_eq!(products[0].read_id, "read");
        assert_eq!(products[0].output_header, "@read/1 extra metadata");
        assert_eq!(products[0].merged_len, 4);
        assert_eq!(products[0].min_qual, 0);
        assert_eq!(products[0].max_qual, 40);
        assert!((products[0].avg_qual - 17.5).abs() < f64::EPSILON);
        assert_eq!(
            products[0].sequence_hash,
            "sha256:1dff3e84fe7877e0673b69bbddcf40124e396e3f9943dd890c91b6a09adb9af0"
        );
        assert_eq!(products[1].read_id, "other/10");
    }

    #[test]
    fn reads_gzip_fastq_products() {
        let path = write_temp_gzip_fastq("@read/2\nAC\n+\nII\n");

        let products = read_merged_products(&path).expect("gzip FASTQ should parse");

        assert_eq!(products.len(), 1);
        assert_eq!(products[0].read_id, "read");
        assert_eq!(products[0].merged_len, 2);
    }

    #[test]
    fn missing_merged_fastq_is_an_error() {
        let path =
            std::env::temp_dir().join(format!("pairasm-products-missing-{}.fastq", Uuid::new_v4()));

        let error = match read_merged_products(&path) {
            Ok(_) => panic!("missing FASTQ should fail"),
            Err(error) => error,
        };

        assert!(error.to_string().contains("merged FASTQ does not exist"));
    }

    #[test]
    fn rejects_malformed_fastq_products() {
        for (name, contents) in [
            ("bad-header", "read\nA\n+\n!\n"),
            ("bad-separator", "@read\nA\n-\n!\n"),
            ("length-mismatch", "@read\nAC\n+\n!\n"),
            ("invalid-quality", "@read\nA\n+\n \n"),
            ("truncated", "@read\nA\n+\n"),
            ("duplicate-read-id", "@read/1\nA\n+\n!\n@read/2\nA\n+\n!\n"),
        ] {
            let path = write_named_temp_fastq(name, contents);

            let result = read_merged_products(&path);

            assert!(result.is_err(), "{name} should fail");
        }
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

    fn write_temp_fastq(contents: &str) -> PathBuf {
        write_named_temp_fastq("valid", contents)
    }

    fn write_named_temp_fastq(name: &str, contents: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("pairasm-products-test-{}", Uuid::new_v4()));
        fs::create_dir(&dir).expect("temp test directory should be created");
        let path = dir.join(format!("{name}.fastq"));
        fs::write(&path, contents).expect("temp FASTQ should be written");
        path
    }

    fn write_temp_gzip_fastq(contents: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("pairasm-products-test-{}", Uuid::new_v4()));
        fs::create_dir(&dir).expect("temp test directory should be created");
        let path = dir.join("valid.fastq.gz");
        let file = fs::File::create(&path).expect("temp gzip FASTQ should be created");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder
            .write_all(contents.as_bytes())
            .expect("gzip FASTQ contents should be written");
        encoder.finish().expect("gzip FASTQ should finish");
        path
    }
}
