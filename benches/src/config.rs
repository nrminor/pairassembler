use std::{
    collections::BTreeMap,
    env,
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use color_eyre::eyre::{Result, WrapErr, bail};

use crate::model::{Dataset, SourceMetadata, SubsetMetadata, Tool, ToolPaths};

pub fn check_tools() -> Result<()> {
    require_command("curl")?;
    let paths = ToolPaths::from_environment()?;
    let tools = [
        ("pairasm", paths.pairasm.as_path()),
        ("fastp", paths.fastp.as_path()),
        ("bbmerge", paths.bbmerge.as_path()),
        ("vsearch", paths.vsearch.as_path()),
        ("hyperfine", paths.hyperfine.as_path()),
    ];

    for (name, path) in tools {
        print_version(name, path)?;
    }

    Ok(())
}

impl ToolPaths {
    pub fn from_environment() -> Result<Self> {
        let file_env = read_env_files()?;
        Ok(Self {
            pairasm: resolve_pairasm(&file_env)?,
            fastp: resolve_tool("FASTP_BIN", "fastp", &file_env)?,
            bbmerge: resolve_tool_any("BBMERGE_BIN", &["bbmerge", "bbmerge.sh"], &file_env)?,
            vsearch: resolve_tool("VSEARCH_BIN", "vsearch", &file_env)?,
            hyperfine: resolve_tool("HYPERFINE_BIN", "hyperfine", &file_env)?,
        })
    }

    pub fn path_for(&self, tool: Tool) -> &Path {
        match tool {
            Tool::Pairasm => &self.pairasm,
            Tool::Fastp => &self.fastp,
            Tool::Bbmerge => &self.bbmerge,
            Tool::Vsearch => &self.vsearch,
        }
    }
}

pub fn read_datasets(path: &Path) -> Result<Vec<Dataset>> {
    let file = File::open(path).wrap_err_with(|| format!("failed to open {}", path.display()))?;
    let mut datasets = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split('\t');
        let name = fields.next().unwrap_or_default().to_owned();
        let accession = fields.next().unwrap_or_default().to_owned();
        if name.is_empty() || accession.is_empty() {
            bail!("dataset rows must have at least name and accession: {trimmed}");
        }
        let default_read_pairs = fields.next().and_then(|raw| raw.parse().ok());
        let note = fields.next().unwrap_or_default().to_owned();
        datasets.push(Dataset {
            name,
            accession,
            default_read_pairs,
            note,
        });
    }
    Ok(datasets)
}

pub fn read_source_metadata(data_root: &Path, name: &str) -> Result<SourceMetadata> {
    let path = data_root.join("raw").join(name).join("source.tsv");
    let row = first_tsv_data_row(&path)?;
    if row.len() < 10 {
        bail!(
            "source metadata row had too few columns: {}",
            path.display()
        );
    }
    Ok(SourceMetadata {
        name: row[0].clone(),
        accession: row[1].clone(),
        r1: PathBuf::from(&row[8]),
        r2: PathBuf::from(&row[9]),
    })
}

pub fn read_subset_metadata(
    data_root: &Path,
    name: &str,
    read_pairs: usize,
) -> Result<SubsetMetadata> {
    let path = data_root
        .join("subset")
        .join(name)
        .join(format!("{read_pairs}_pairs"))
        .join("subset.tsv");
    let row = first_tsv_data_row(&path)?;
    if row.len() < 5 {
        bail!(
            "subset metadata row had too few columns; rerun bench-real-prepare: {}",
            path.display()
        );
    }
    Ok(SubsetMetadata {
        name: row[0].clone(),
        accession: row[1].clone(),
        read_pairs: row[2].parse()?,
        r1: PathBuf::from(&row[3]),
        r2: PathBuf::from(&row[4]),
    })
}

pub fn first_tsv_data_row(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path)?;
    for (index, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if index == 0 || line.trim().is_empty() {
            continue;
        }
        return Ok(line.split('\t').map(str::to_owned).collect());
    }
    bail!("no data rows found in {}", path.display())
}

pub fn write_tool_versions(run_dir: &Path, paths: &ToolPaths) -> Result<()> {
    let mut writer = BufWriter::new(File::create(
        run_dir.join("metadata").join("tool_versions.tsv"),
    )?);
    writeln!(writer, "tool\tpath\tversion")?;
    for (name, path) in [
        ("pairasm", paths.pairasm.as_path()),
        ("fastp", paths.fastp.as_path()),
        ("bbmerge", paths.bbmerge.as_path()),
        ("vsearch", paths.vsearch.as_path()),
        ("hyperfine", paths.hyperfine.as_path()),
    ] {
        let version =
            version_string(path).unwrap_or_else(|error| format!("version unavailable: {error}"));
        writeln!(writer, "{name}\t{}\t{version}", path.display())?;
    }
    writer.flush()?;
    Ok(())
}

pub fn require_command(binary: &str) -> Result<()> {
    find_on_path(OsStr::new(binary))
        .map(|_| ())
        .ok_or_else(|| color_eyre::eyre::eyre!("required command not found on PATH: {binary}"))
}

fn read_env_files() -> Result<BTreeMap<String, String>> {
    let mut values = BTreeMap::new();
    for path in ["benches/config/benchmark.env", "benches/config/tools.env"] {
        let path = Path::new(path);
        if path.exists() {
            read_env_file(path, &mut values)?;
        }
    }
    Ok(values)
}

fn read_env_file(path: &Path, values: &mut BTreeMap<String, String>) -> Result<()> {
    let file = File::open(path)?;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let Some((key, value)) = trimmed.split_once('=') else {
            continue;
        };
        values.insert(key.trim().to_owned(), unquote(value.trim()).to_owned());
    }
    Ok(())
}

fn unquote(value: &str) -> &str {
    value
        .strip_prefix('"')
        .and_then(|without_prefix| without_prefix.strip_suffix('"'))
        .unwrap_or(value)
}

fn env_or_file(key: &str, file_env: &BTreeMap<String, String>) -> Option<String> {
    env::var(key).ok().or_else(|| file_env.get(key).cloned())
}

fn resolve_pairasm(file_env: &BTreeMap<String, String>) -> Result<PathBuf> {
    if let Some(path) = env_or_file("PAIRASM_BIN", file_env) {
        let path = PathBuf::from(path);
        if path.exists() || find_on_path(path.as_os_str()).is_some() {
            return Ok(path);
        }
    }

    let release_path = PathBuf::from("target/release/pairasm");
    if release_path.exists() {
        return Ok(release_path);
    }

    resolve_tool("PAIRASM_BIN", "pairasm", file_env)
}

fn resolve_tool(
    key: &str,
    default_name: &str,
    file_env: &BTreeMap<String, String>,
) -> Result<PathBuf> {
    let configured = env_or_file(key, file_env).unwrap_or_else(|| default_name.to_owned());
    let path = PathBuf::from(&configured);
    if path.components().count() > 1 && path.exists() {
        return Ok(path);
    }
    find_on_path(path.as_os_str())
        .or_else(|| path.exists().then_some(path.clone()))
        .ok_or_else(|| color_eyre::eyre::eyre!("required benchmark tool not found: {configured}"))
}

fn resolve_tool_any(
    key: &str,
    default_names: &[&str],
    file_env: &BTreeMap<String, String>,
) -> Result<PathBuf> {
    if let Some(configured) = env_or_file(key, file_env) {
        return resolve_tool(key, &configured, file_env);
    }

    for default_name in default_names {
        if let Some(path) = find_on_path(OsStr::new(default_name)) {
            return Ok(path);
        }
    }

    bail!(
        "required benchmark tool not found: set {key} or install one of: {}",
        default_names.join(", ")
    )
}

fn find_on_path(binary: &OsStr) -> Option<PathBuf> {
    let path_var = env::var_os("PATH")?;
    env::split_paths(&path_var)
        .map(|dir| dir.join(binary))
        .find(|candidate| candidate.exists())
}

fn print_version(name: &str, path: &Path) -> Result<()> {
    let output = Command::new(path)
        .arg("--version")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .wrap_err_with(|| format!("failed to run {} --version", path.display()))?;
    let mut version = String::from_utf8_lossy(&output.stdout).trim().to_owned();
    if version.is_empty() {
        version = String::from_utf8_lossy(&output.stderr).trim().to_owned();
    }
    println!("{name}\t{}\t{version}", path.display());
    Ok(())
}

fn version_string(path: &Path) -> Result<String> {
    let output = Command::new(path).arg("--version").output()?;
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_owned();
    if !stdout.is_empty() {
        return Ok(stdout);
    }
    Ok(String::from_utf8_lossy(&output.stderr).trim().to_owned())
}
