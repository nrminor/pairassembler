use std::{fs, path::Path, process::Command};

use color_eyre::eyre::{Result, bail};

use crate::{
    cli::CommonOptions,
    config::{Dataset, first_tsv_data_row, read_datasets, require_command},
    fastq::validate_gzip,
    process::run_command,
    ui,
};

pub(crate) fn fetch_ena(options: &CommonOptions) -> Result<()> {
    require_command("curl")?;
    let datasets = read_datasets(&options.config)?;
    let raw_root = options.data_root.join("raw");
    fs::create_dir_all(&raw_root)?;

    for dataset in datasets {
        fetch_dataset(&dataset, &raw_root)?;
    }

    Ok(())
}

fn fetch_dataset(dataset: &Dataset, raw_root: &Path) -> Result<()> {
    let out_dir = raw_root.join(&dataset.name);
    fs::create_dir_all(&out_dir)?;
    let report_path = out_dir.join("ena_filereport.tsv");
    let report_url = format!(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=run_accession,instrument_platform,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_md5&format=tsv&download=true",
        dataset.accession
    );

    eprintln!(
        "{} {} {}",
        ui::muted_stderr("Resolving ENA FASTQs for"),
        ui::dataset_stderr(&dataset.name),
        ui::muted_stderr(format!("({})", dataset.accession))
    );
    run_command(
        Command::new("curl")
            .args(["-fsSL", &report_url, "-o"])
            .arg(&report_path),
    )?;

    let row = first_tsv_data_row(&report_path)?;
    if row.len() < 8 {
        bail!(
            "ENA filereport row for {} had too few columns",
            dataset.accession
        );
    }
    let layout = &row[2];
    if layout != "PAIRED" {
        bail!(
            "{} is not paired according to ENA: {layout}",
            dataset.accession
        );
    }

    let urls: Vec<&str> = row[6].split(';').collect();
    if urls.len() < 2 {
        bail!(
            "ENA did not report two FASTQ URLs for {}",
            dataset.accession
        );
    }
    let r1_path = out_dir.join(format!("{}_1.fastq.gz", dataset.accession));
    let r2_path = out_dir.join(format!("{}_2.fastq.gz", dataset.accession));
    let r1_url = format!("https://{}", urls[0]);
    let r2_url = format!("https://{}", urls[1]);
    download_if_missing(&r1_url, &r1_path)?;
    download_if_missing(&r2_url, &r2_path)?;
    validate_gzip(&r1_path)?;
    validate_gzip(&r2_path)?;

    let source_path = out_dir.join("source.tsv");
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(source_path)?;
    writer.write_record([
        "name",
        "accession",
        "run",
        "platform",
        "layout",
        "strategy",
        "read_count",
        "base_count",
        "r1",
        "r2",
        "note",
    ])?;
    writer.write_record([
        dataset.name.clone(),
        dataset.accession.clone(),
        row[0].clone(),
        row[1].clone(),
        row[2].clone(),
        row[3].clone(),
        row[4].clone(),
        row[5].clone(),
        r1_path.to_string_lossy().into_owned(),
        r2_path.to_string_lossy().into_owned(),
        dataset.note.clone(),
    ])?;
    writer.flush()?;
    Ok(())
}

fn download_if_missing(url: &str, path: &Path) -> Result<()> {
    if path.exists() && path.metadata()?.len() > 0 {
        eprintln!(
            "{} {}",
            ui::muted_stderr("Using cached"),
            ui::path_stderr(path.display())
        );
        return Ok(());
    }
    eprintln!("{} {url}", ui::muted_stderr("Downloading"));
    run_command(Command::new("curl").args(["-fL", url, "-o"]).arg(path))
}
