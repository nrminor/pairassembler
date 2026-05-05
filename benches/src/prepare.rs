use std::{
    fs::{self, File},
    io::{BufWriter, Write},
};

use color_eyre::eyre::Result;

use crate::{
    cli::PrepareOptions,
    config::{read_datasets, read_source_metadata},
    fastq::{fastq_record_count, file_size, validate_gzip, write_first_fastq_records},
    model::SourceMetadata,
};

pub fn prepare_subsets(options: &PrepareOptions) -> Result<()> {
    let datasets = read_datasets(&options.common.config)?;

    for dataset in datasets {
        let requested_pairs = dataset
            .default_read_pairs
            .unwrap_or(options.read_pairs)
            .min(options.read_pairs);
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
            "Preparing first {read_pairs} read pairs for {}",
            source.name
        );
        write_first_fastq_records(&source.r1, &out_r1, read_pairs)?;
        write_first_fastq_records(&source.r2, &out_r2, read_pairs)?;
    }
    validate_gzip(&out_r1)?;
    validate_gzip(&out_r2)?;

    let mut writer = BufWriter::new(File::create(out_dir.join("subset.tsv"))?);
    writeln!(
        writer,
        "name\taccession\tread_pairs\tr1\tr2\tr1_bytes\tr2_bytes\tr1_records\tr2_records"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        source.name,
        source.accession,
        read_pairs,
        out_r1.display(),
        out_r2.display(),
        file_size(&out_r1)?,
        file_size(&out_r2)?,
        fastq_record_count(&out_r1)?,
        fastq_record_count(&out_r2)?
    )?;
    writer.flush()?;
    Ok(())
}
