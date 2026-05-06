use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use color_eyre::eyre::{Result, bail};

use crate::{cli::SummarizeOptions, compare::write_pairwise_agreement};

pub fn summarize(options: &SummarizeOptions) -> Result<()> {
    let run_dir = match &options.run_dir {
        Some(path) => path.clone(),
        None if options.latest => latest_run_dir(&options.runs_root)?,
        None => bail!("provide --run-dir or use --latest"),
    };

    let result_files = collect_files_named(&run_dir, "result.tsv")?;
    if result_files.is_empty() {
        bail!("no result.tsv files found under {}", run_dir.display());
    }

    let summary_dir = run_dir.join("summary");
    fs::create_dir_all(&summary_dir)?;
    let summary_path = summary_dir.join("results.tsv");
    let mut writer = BufWriter::new(File::create(&summary_path)?);
    let mut wrote_header = false;

    for path in &result_files {
        let file = File::open(path)?;
        let mut lines = BufReader::new(file).lines();
        if let Some(header) = lines.next() {
            let header = header?;
            if !wrote_header {
                writeln!(writer, "{header}")?;
                wrote_header = true;
            }
        }
        for line in lines {
            writeln!(writer, "{}", line?)?;
        }
    }

    writer.flush()?;
    let agreement = write_pairwise_agreement(&result_files, &summary_dir)?;
    println!("{}", summary_path.display());
    println!("{}", agreement.tsv_path.display());
    println!("{}", agreement.markdown_path.display());
    Ok(())
}

fn latest_run_dir(runs_root: &Path) -> Result<PathBuf> {
    let mut dirs = fs::read_dir(runs_root)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.is_dir())
        .collect::<Vec<_>>();
    dirs.sort();
    dirs.pop().ok_or_else(|| {
        color_eyre::eyre::eyre!("no benchmark runs found under {}", runs_root.display())
    })
}

fn collect_files_named(root: &Path, file_name: &str) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    collect_files_named_inner(root, file_name, &mut files)?;
    files.sort();
    Ok(files)
}

fn collect_files_named_inner(root: &Path, file_name: &str, files: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            collect_files_named_inner(&path, file_name, files)?;
        } else if path.file_name().is_some_and(|name| name == file_name) {
            files.push(path);
        }
    }
    Ok(())
}
