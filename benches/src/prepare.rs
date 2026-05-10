use std::fs;

use color_eyre::eyre::Result;

use crate::{
    cli::PrepareOptions,
    config::{effective_read_pairs, read_datasets, read_source_metadata},
    fastq::{fastq_record_count, file_size, validate_gzip, write_first_fastq_records},
    model::SourceMetadata,
    ui,
};

pub fn prepare_subsets(options: &PrepareOptions) -> Result<()> {
    let datasets = read_datasets(&options.common.config)?;

    for dataset in datasets {
        let requested_pairs = effective_read_pairs(&dataset, options.read_pairs);
        let source = read_source_metadata(&options.common.data_root, &dataset.name)?;
        prepare_subset(&source, requested_pairs, &options.common.data_root)?;
    }

    Ok(())
}

fn prepare_subset(
    source: &SourceMetadata,
    read_pairs: usize,
    data_root: &std::path::Path,
) -> Result<()> {
    let out_dir = data_root
        .join("subset")
        .join(&source.name)
        .join(format!("{read_pairs}_pairs"));
    fs::create_dir_all(&out_dir)?;
    let out_r1 = out_dir.join(format!(
        "{}_{}_pairs_1.fastq.gz",
        source.accession, read_pairs
    ));
    let out_r2 = out_dir.join(format!(
        "{}_{}_pairs_2.fastq.gz",
        source.accession, read_pairs
    ));
    if !out_r1.exists() || !out_r2.exists() {
        eprintln!(
            "{} {} {}",
            ui::muted_stderr("Preparing first"),
            ui::muted_stderr(format!("{read_pairs} read pairs for")),
            ui::dataset_stderr(&source.name)
        );
        write_first_fastq_records(&source.r1, &out_r1, read_pairs)?;
        write_first_fastq_records(&source.r2, &out_r2, read_pairs)?;
    } else {
        eprintln!(
            "{} {} {}",
            ui::muted_stderr("Using cached"),
            ui::muted_stderr(format!("{read_pairs}-pair subset for")),
            ui::dataset_stderr(&source.name)
        );
    }
    validate_gzip(&out_r1)?;
    validate_gzip(&out_r2)?;

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out_dir.join("subset.tsv"))?;
    writer.write_record([
        "name",
        "accession",
        "read_pairs",
        "r1",
        "r2",
        "r1_bytes",
        "r2_bytes",
        "r1_records",
        "r2_records",
    ])?;
    writer.write_record([
        source.name.clone(),
        source.accession.clone(),
        read_pairs.to_string(),
        out_r1.to_string_lossy().into_owned(),
        out_r2.to_string_lossy().into_owned(),
        file_size(&out_r1)?.to_string(),
        file_size(&out_r2)?.to_string(),
        fastq_record_count(&out_r1)?.to_string(),
        fastq_record_count(&out_r2)?.to_string(),
    ])?;
    writer.flush()?;
    Ok(())
}
