use std::{fs, path::Path};

use color_eyre::eyre::{Result, WrapErr, bail};

use crate::{config::SubsetMetadata, fastq::fastq_record_count, tool::Tool};

pub(crate) fn validate_tool_run(
    tool: Tool,
    _mode: crate::cli::BenchmarkMode,
    subset: &SubsetMetadata,
    out_dir: &Path,
    merged_output: &Path,
    merged_reads: usize,
) -> Result<()> {
    match tool {
        Tool::Bbmerge => validate_bbmerge_run(subset, out_dir, merged_output, merged_reads),
        Tool::Vsearch => validate_vsearch_run(subset, out_dir, merged_output, merged_reads),
        Tool::Pairasm | Tool::Fastp => Ok(()),
    }
}

fn validate_bbmerge_run(
    subset: &SubsetMetadata,
    out_dir: &Path,
    merged_output: &Path,
    merged_reads: usize,
) -> Result<()> {
    let stderr = read_tool_log(out_dir, "bbmerge", "stderr")?;
    reject_java_exceptions("bbmerge", out_dir, &stderr)?;

    let reported_pairs = parse_colon_metric(&stderr, "Pairs")?;
    let joined = parse_colon_metric(&stderr, "Joined")?;
    let no_solution = parse_colon_metric(&stderr, "No Solution")?;
    let too_short = parse_colon_metric(&stderr, "Too Short")?;

    if reported_pairs != subset.read_pairs {
        bail!(
            "bbmerge processed {reported_pairs} pairs, expected {}; see {}",
            subset.read_pairs,
            out_dir.display()
        );
    }
    if joined != merged_reads {
        bail!(
            "bbmerge reported {joined} joined pairs, but {} contains {merged_reads} records",
            merged_output.display()
        );
    }
    if joined + no_solution + too_short != subset.read_pairs {
        bail!(
            "bbmerge accounting does not sum to expected pairs: joined={joined}, no_solution={no_solution}, too_short={too_short}, expected={}",
            subset.read_pairs
        );
    }
    validate_bbmerge_unmerged_output(subset, out_dir, joined)
}

fn validate_bbmerge_unmerged_output(
    subset: &SubsetMetadata,
    out_dir: &Path,
    joined: usize,
) -> Result<()> {
    let expected_pairs = subset.read_pairs - joined;
    for name in ["bbmerge.unmerged1.fastq", "bbmerge.unmerged2.fastq"] {
        let path = out_dir.join(name);
        if path.exists() {
            let records = fastq_record_count(&path)?;
            if records != expected_pairs {
                bail!(
                    "{} contains {records} records, expected {expected_pairs}",
                    path.display()
                );
            }
        }
    }
    Ok(())
}

fn validate_vsearch_run(
    subset: &SubsetMetadata,
    out_dir: &Path,
    merged_output: &Path,
    merged_reads: usize,
) -> Result<()> {
    let stderr = read_tool_log(out_dir, "vsearch", "stderr")?;
    let reported_pairs = parse_vsearch_metric(&stderr, &["Pairs"])?;
    let merged = parse_vsearch_metric(&stderr, &["Merged"])?;
    let not_merged = parse_vsearch_metric(&stderr, &["Not", "merged"])?;

    if reported_pairs != subset.read_pairs {
        bail!(
            "vsearch processed {reported_pairs} pairs, expected {}; see {}",
            subset.read_pairs,
            out_dir.display()
        );
    }
    if merged != merged_reads {
        bail!(
            "vsearch reported {merged} merged pairs, but {} contains {merged_reads} records",
            merged_output.display()
        );
    }
    if merged + not_merged != subset.read_pairs {
        bail!(
            "vsearch accounting does not sum to expected pairs: merged={merged}, not_merged={not_merged}, expected={}",
            subset.read_pairs
        );
    }

    Ok(())
}

fn read_tool_log(out_dir: &Path, tool_name: &str, stream_name: &str) -> Result<String> {
    let path = out_dir.join(format!("{tool_name}.{stream_name}.log"));
    fs::read_to_string(&path).wrap_err_with(|| format!("failed to read {}", path.display()))
}

fn reject_java_exceptions(tool_name: &str, out_dir: &Path, stderr: &str) -> Result<()> {
    if stderr.contains("Exception in thread") || stderr.contains("java.lang.") {
        bail!(
            "{tool_name} emitted a Java exception; refusing to summarize invalid benchmark artifact under {}",
            out_dir.display()
        );
    }
    Ok(())
}

fn parse_colon_metric(text: &str, label: &str) -> Result<usize> {
    let prefix = format!("{label}:");
    for line in text.lines() {
        let trimmed = line.trim_start();
        if let Some(rest) = trimmed.strip_prefix(&prefix) {
            let raw = rest
                .split_whitespace()
                .next()
                .ok_or_else(|| color_eyre::eyre::eyre!("missing value for {label}"))?;
            return raw.parse::<usize>().wrap_err_with(|| {
                format!("failed to parse {label} metric value {raw:?} as an integer")
            });
        }
    }
    bail!("missing {label} metric")
}

fn parse_vsearch_metric(text: &str, label_words: &[&str]) -> Result<usize> {
    for line in text.lines() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        let label_matches = fields
            .iter()
            .skip(1)
            .take(label_words.len())
            .copied()
            .eq(label_words.iter().copied());
        if fields.len() > label_words.len() && label_matches {
            return fields[0].parse::<usize>().wrap_err_with(|| {
                format!(
                    "failed to parse vsearch {} metric value {:?} as an integer",
                    label_words.join(" "),
                    fields[0]
                )
            });
        }
    }
    bail!("missing vsearch {} metric", label_words.join(" "))
}

#[cfg(test)]
mod tests {
    use super::{parse_colon_metric, parse_vsearch_metric};

    #[test]
    fn parses_bbmerge_colon_metrics() {
        let stderr = "\
            Pairs:             100\n\
            Joined:            82\n\
            No Solution:       17\n\
            Too Short:         1\n";

        assert_eq!(parse_colon_metric(stderr, "Pairs").unwrap(), 100);
        assert_eq!(parse_colon_metric(stderr, "Joined").unwrap(), 82);
        assert_eq!(parse_colon_metric(stderr, "No Solution").unwrap(), 17);
        assert_eq!(parse_colon_metric(stderr, "Too Short").unwrap(), 1);
    }

    #[test]
    fn rejects_missing_or_invalid_bbmerge_colon_metrics() {
        assert!(parse_colon_metric("Joined: 10\n", "Pairs").is_err());
        assert!(parse_colon_metric("Pairs: nope\n", "Pairs").is_err());
    }

    #[test]
    fn parses_vsearch_metrics() {
        let stderr = "\
            100  Pairs\n\
             82  Merged (82.0%)\n\
             18  Not merged (18.0%)\n";

        assert_eq!(parse_vsearch_metric(stderr, &["Pairs"]).unwrap(), 100);
        assert_eq!(parse_vsearch_metric(stderr, &["Merged"]).unwrap(), 82);
        assert_eq!(
            parse_vsearch_metric(stderr, &["Not", "merged"]).unwrap(),
            18
        );
    }

    #[test]
    fn rejects_missing_or_invalid_vsearch_metrics() {
        assert!(parse_vsearch_metric("82 Merged\n", &["Pairs"]).is_err());
        assert!(parse_vsearch_metric("nope Pairs\n", &["Pairs"]).is_err());
    }
}
