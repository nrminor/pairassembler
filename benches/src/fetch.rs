use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    process::Command,
};

use color_eyre::eyre::{Result, bail};

use crate::{
    cli::CommonOptions,
    config::{first_tsv_data_row, read_datasets, require_command},
    fastq::validate_gzip,
    model::Dataset,
    process::run_command,
    ui,
};

pub fn fetch_ena(options: &CommonOptions) -> Result<()> {
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
    download_if_missing(&format!("https://{}", urls[0]), &r1_path)?;
    download_if_missing(&format!("https://{}", urls[1]), &r2_path)?;
    validate_gzip(&r1_path)?;
    validate_gzip(&r2_path)?;

    let source_path = out_dir.join("source.tsv");
    let mut writer = BufWriter::new(File::create(source_path)?);
    writeln!(
        writer,
        "name\taccession\trun\tplatform\tlayout\tstrategy\tread_count\tbase_count\tr1\tr2\tnote"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        dataset.name,
        dataset.accession,
        row[0],
        row[1],
        row[2],
        row[3],
        row[4],
        row[5],
        r1_path.display(),
        r2_path.display(),
        dataset.note
    )?;
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
