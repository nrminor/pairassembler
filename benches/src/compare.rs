use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use color_eyre::eyre::{Context, Result, bail, eyre};

pub struct AgreementArtifacts {
    pub tsv_path: PathBuf,
    pub markdown_path: PathBuf,
}

pub fn write_pairwise_agreement(
    result_files: &[PathBuf],
    summary_dir: &Path,
) -> Result<AgreementArtifacts> {
    let tool_results = read_tool_results(result_files)?;
    let manifests = read_manifests(&tool_results)?;
    let agreements = compute_pairwise_agreements(&tool_results, &manifests);

    let tsv_path = summary_dir.join("pairwise-agreement.tsv");
    write_pairwise_agreement_tsv(&tsv_path, &agreements)?;

    let markdown_path = summary_dir.join("pairwise-agreement.md");
    write_pairwise_agreement_markdown(&markdown_path, &agreements)?;

    Ok(AgreementArtifacts {
        tsv_path,
        markdown_path,
    })
}

fn read_tool_results(result_files: &[PathBuf]) -> Result<Vec<ToolResult>> {
    let mut results = Vec::new();
    for path in result_files {
        let file =
            File::open(path).wrap_err_with(|| format!("failed to open {}", path.display()))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| eyre!("{} is empty", path.display()))??;
        let columns = HeaderColumns::from_header(&header)?;

        for line in lines {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            results.push(ToolResult {
                dataset: columns.field(&fields, "dataset")?.to_owned(),
                tool: columns.field(&fields, "tool")?.to_owned(),
                merged_reads: parse_usize(columns.field(&fields, "merged_reads")?, "merged_reads")?,
                manifest_rows: parse_usize(
                    columns.field(&fields, "manifest_rows")?,
                    "manifest_rows",
                )?,
                output_dir: PathBuf::from(columns.field(&fields, "output_dir")?),
            });
        }
    }
    results.sort_by(|a, b| a.dataset.cmp(&b.dataset).then_with(|| a.tool.cmp(&b.tool)));
    Ok(results)
}

fn read_manifests(
    tool_results: &[ToolResult],
) -> Result<BTreeMap<(String, String), BTreeMap<String, ManifestRecord>>> {
    let mut manifests = BTreeMap::new();
    for result in tool_results {
        let path = result.output_dir.join("merged-manifest.tsv");
        let records = read_manifest(&path)?;
        if records.len() != result.manifest_rows {
            bail!(
                "{} contained {} rows, but result.tsv reported {} manifest rows",
                path.display(),
                records.len(),
                result.manifest_rows
            );
        }
        manifests.insert((result.dataset.clone(), result.tool.clone()), records);
    }
    Ok(manifests)
}

fn read_manifest(path: &Path) -> Result<BTreeMap<String, ManifestRecord>> {
    let file = File::open(path).wrap_err_with(|| format!("failed to open {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| eyre!("{} is empty", path.display()))??;
    let columns = HeaderColumns::from_header(&header)?;

    let mut records = BTreeMap::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let read_id = columns.field(&fields, "read_id")?.to_owned();
        let record = ManifestRecord {
            merged_len: parse_usize(columns.field(&fields, "merged_len")?, "merged_len")?,
            sequence_hash: columns.field(&fields, "sequence_hash")?.to_owned(),
            quality_hash: columns.field(&fields, "quality_hash")?.to_owned(),
        };
        if records.insert(read_id.clone(), record).is_some() {
            bail!("duplicate read_id {read_id:?} in {}", path.display());
        }
    }
    Ok(records)
}

fn compute_pairwise_agreements(
    tool_results: &[ToolResult],
    manifests: &BTreeMap<(String, String), BTreeMap<String, ManifestRecord>>,
) -> Vec<PairwiseAgreement> {
    let mut by_dataset: BTreeMap<&str, Vec<&ToolResult>> = BTreeMap::new();
    for result in tool_results {
        by_dataset.entry(&result.dataset).or_default().push(result);
    }

    let mut agreements = Vec::new();
    for (dataset, results) in by_dataset {
        for left_idx in 0..results.len() {
            for right_idx in (left_idx + 1)..results.len() {
                let left = results[left_idx];
                let right = results[right_idx];
                let left_records = &manifests[&(left.dataset.clone(), left.tool.clone())];
                let right_records = &manifests[&(right.dataset.clone(), right.tool.clone())];
                agreements.push(compare_pair(
                    dataset,
                    left,
                    right,
                    left_records,
                    right_records,
                ));
            }
        }
    }
    agreements
}

fn compare_pair(
    dataset: &str,
    left: &ToolResult,
    right: &ToolResult,
    left_records: &BTreeMap<String, ManifestRecord>,
    right_records: &BTreeMap<String, ManifestRecord>,
) -> PairwiseAgreement {
    let left_ids = left_records.keys().cloned().collect::<BTreeSet<_>>();
    let right_ids = right_records.keys().cloned().collect::<BTreeSet<_>>();
    let merged_both_ids = left_ids.intersection(&right_ids).collect::<Vec<_>>();
    let only_left = left_ids.difference(&right_ids).count();
    let only_right = right_ids.difference(&left_ids).count();
    let union = left_ids.union(&right_ids).count();

    let mut same_sequence = 0usize;
    let mut different_sequence = 0usize;
    let mut same_quality = 0usize;
    let mut different_quality = 0usize;
    let mut same_length_different_sequence = 0usize;
    let mut different_length = 0usize;

    for read_id in &merged_both_ids {
        let left_record = &left_records[*read_id];
        let right_record = &right_records[*read_id];
        if left_record.sequence_hash == right_record.sequence_hash {
            same_sequence += 1;
        } else {
            different_sequence += 1;
            if left_record.merged_len == right_record.merged_len {
                same_length_different_sequence += 1;
            } else {
                different_length += 1;
            }
        }

        if left_record.quality_hash == right_record.quality_hash {
            same_quality += 1;
        } else {
            different_quality += 1;
        }
    }

    PairwiseAgreement {
        dataset: dataset.to_owned(),
        tool_a: left.tool.clone(),
        tool_b: right.tool.clone(),
        merged_a: left.merged_reads,
        merged_b: right.merged_reads,
        merged_both: merged_both_ids.len(),
        only_a: only_left,
        only_b: only_right,
        jaccard_merged: jaccard(merged_both_ids.len(), union),
        same_sequence,
        different_sequence,
        same_quality,
        different_quality,
        same_length_different_sequence,
        different_length,
    }
}

fn write_pairwise_agreement_tsv(path: &Path, agreements: &[PairwiseAgreement]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "dataset\ttool_a\ttool_b\tmerged_a\tmerged_b\tmerged_both\tonly_a\tonly_b\tjaccard_merged\tsame_sequence\tdifferent_sequence\tsame_quality\tdifferent_quality\tsame_length_different_sequence\tdifferent_length"
    )?;
    for agreement in agreements {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
            agreement.dataset,
            agreement.tool_a,
            agreement.tool_b,
            agreement.merged_a,
            agreement.merged_b,
            agreement.merged_both,
            agreement.only_a,
            agreement.only_b,
            agreement.jaccard_merged,
            agreement.same_sequence,
            agreement.different_sequence,
            agreement.same_quality,
            agreement.different_quality,
            agreement.same_length_different_sequence,
            agreement.different_length
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_pairwise_agreement_markdown(path: &Path, agreements: &[PairwiseAgreement]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# Pairwise merged-read agreement")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Agreement is computed from merged-read manifests. It compares which read IDs were merged and whether commonly merged products have identical sequence or quality hashes."
    )?;

    let mut current_dataset = None::<&str>;
    for agreement in agreements {
        if current_dataset != Some(&agreement.dataset) {
            current_dataset = Some(&agreement.dataset);
            writeln!(writer)?;
            writeln!(writer, "## {}", agreement.dataset)?;
            writeln!(writer)?;
            writeln!(
                writer,
                "| tool A | tool B | merged A | merged B | both | only A | only B | Jaccard | same seq | diff seq | same qual | diff qual |"
            )?;
            writeln!(
                writer,
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
            )?;
        }
        writeln!(
            writer,
            "| {} | {} | {} | {} | {} | {} | {} | {:.3} | {} | {} | {} | {} |",
            agreement.tool_a,
            agreement.tool_b,
            agreement.merged_a,
            agreement.merged_b,
            agreement.merged_both,
            agreement.only_a,
            agreement.only_b,
            agreement.jaccard_merged,
            agreement.same_sequence,
            agreement.different_sequence,
            agreement.same_quality,
            agreement.different_quality
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn parse_usize(raw: &str, label: &str) -> Result<usize> {
    raw.parse::<usize>()
        .wrap_err_with(|| format!("failed to parse {label} value {raw:?} as usize"))
}

fn jaccard(intersection: usize, union: usize) -> f64 {
    if union == 0 {
        return 1.0;
    }
    intersection as f64 / union as f64
}

struct HeaderColumns {
    indices: BTreeMap<String, usize>,
}

impl HeaderColumns {
    fn from_header(header: &str) -> Result<Self> {
        let indices = header
            .split('\t')
            .enumerate()
            .map(|(idx, column)| (column.to_owned(), idx))
            .collect::<BTreeMap<_, _>>();
        Ok(Self { indices })
    }

    fn field<'a>(&self, fields: &'a [&str], column: &str) -> Result<&'a str> {
        let index = self
            .indices
            .get(column)
            .ok_or_else(|| eyre!("missing required column {column:?}"))?;
        fields
            .get(*index)
            .copied()
            .ok_or_else(|| eyre!("row is missing field for required column {column:?}"))
    }
}

struct ToolResult {
    dataset: String,
    tool: String,
    merged_reads: usize,
    manifest_rows: usize,
    output_dir: PathBuf,
}

struct ManifestRecord {
    merged_len: usize,
    sequence_hash: String,
    quality_hash: String,
}

struct PairwiseAgreement {
    dataset: String,
    tool_a: String,
    tool_b: String,
    merged_a: usize,
    merged_b: usize,
    merged_both: usize,
    only_a: usize,
    only_b: usize,
    jaccard_merged: f64,
    same_sequence: usize,
    different_sequence: usize,
    same_quality: usize,
    different_quality: usize,
    same_length_different_sequence: usize,
    different_length: usize,
}

#[cfg(test)]
mod tests {
    use super::jaccard;

    #[test]
    fn jaccard_is_one_for_empty_sets() {
        assert!((jaccard(0, 0) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn jaccard_is_intersection_over_union() {
        assert!((jaccard(2, 5) - 0.4).abs() < f64::EPSILON);
    }
}
