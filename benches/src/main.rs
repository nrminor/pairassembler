use std::{
    collections::BTreeMap,
    env,
    ffi::OsStr,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use clap::{Parser, Subcommand, ValueEnum};
use color_eyre::eyre::{Result, WrapErr, bail};
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use serde::Deserialize;

const DEFAULT_CONFIG: &str = "benches/config/datasets.tsv";
const DEFAULT_DATA_ROOT: &str = "benches/data";
const DEFAULT_RUNS_ROOT: &str = "benches/runs";
const DEFAULT_READ_PAIRS: usize = 100_000;
const DEFAULT_REPLICATES: usize = 3;
const DEFAULT_THREADS: usize = 8;

#[derive(Debug, Parser)]
#[command(about = "Real-data comparative benchmarks for pairasm")]
struct Cli {
    #[command(subcommand)]
    command: BenchCommand,
}

#[derive(Debug, Subcommand)]
enum BenchCommand {
    /// Check external benchmark tools and print versions.
    Check,
    /// Fetch configured paired FASTQs from ENA.
    Fetch(CommonOptions),
    /// Prepare deterministic first-N-pair FASTQ subsets.
    Prepare(PrepareOptions),
    /// Run pairasm and competitor merge tools through hyperfine.
    Run(RunOptions),
    /// Summarize hyperfine run artifacts into a TSV table.
    Summarize(SummarizeOptions),
}

#[derive(Debug, Clone, Parser)]
struct CommonOptions {
    #[arg(long, default_value = DEFAULT_CONFIG)]
    config: PathBuf,
    #[arg(long, default_value = DEFAULT_DATA_ROOT)]
    data_root: PathBuf,
}

#[derive(Debug, Clone, Parser)]
struct PrepareOptions {
    #[command(flatten)]
    common: CommonOptions,
    #[arg(long, default_value_t = DEFAULT_READ_PAIRS)]
    read_pairs: usize,
}

#[derive(Debug, Clone, Parser)]
struct RunOptions {
    #[command(flatten)]
    common: CommonOptions,
    #[arg(long, default_value = DEFAULT_RUNS_ROOT)]
    runs_root: PathBuf,
    #[arg(long, default_value_t = DEFAULT_READ_PAIRS)]
    read_pairs: usize,
    #[arg(long, default_value_t = DEFAULT_REPLICATES)]
    replicates: usize,
    #[arg(long, default_value_t = DEFAULT_THREADS)]
    threads: usize,
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "pairasm,fastp,bbmerge,vsearch"
    )]
    tools: Vec<Tool>,
    #[arg(long, default_value_t = OutputCompression::Plain)]
    output_compression: OutputCompression,
}

#[derive(Debug, Clone, Parser)]
struct SummarizeOptions {
    #[arg(long, default_value = DEFAULT_RUNS_ROOT)]
    runs_root: PathBuf,
    #[arg(long)]
    run_dir: Option<PathBuf>,
    #[arg(long, default_value_t = true)]
    latest: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Tool {
    Pairasm,
    Fastp,
    Bbmerge,
    Vsearch,
}

impl Tool {
    fn name(self) -> &'static str {
        match self {
            Tool::Pairasm => "pairasm",
            Tool::Fastp => "fastp",
            Tool::Bbmerge => "bbmerge",
            Tool::Vsearch => "vsearch",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum OutputCompression {
    Plain,
    Gzip,
}

impl std::fmt::Display for OutputCompression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Plain => f.write_str("plain"),
            Self::Gzip => f.write_str("gzip"),
        }
    }
}

#[derive(Clone, Debug)]
struct Dataset {
    name: String,
    accession: String,
    default_read_pairs: Option<usize>,
    note: String,
}

#[derive(Clone, Debug)]
struct ToolPaths {
    pairasm: PathBuf,
    fastp: PathBuf,
    bbmerge: PathBuf,
    vsearch: PathBuf,
    hyperfine: PathBuf,
}

#[derive(Debug)]
struct SourceMetadata {
    name: String,
    accession: String,
    r1: PathBuf,
    r2: PathBuf,
}

#[derive(Debug)]
struct SubsetMetadata {
    name: String,
    accession: String,
    read_pairs: usize,
    r1: PathBuf,
    r2: PathBuf,
}

#[derive(Debug)]
struct ToolCommand {
    tool: Tool,
    args: Vec<String>,
    merged_output: PathBuf,
}

#[derive(Debug, Deserialize)]
struct HyperfineReport {
    results: Vec<HyperfineResult>,
}

#[derive(Debug, Deserialize)]
struct HyperfineResult {
    mean: f64,
    stddev: Option<f64>,
    median: f64,
    min: f64,
    max: f64,
    user: f64,
    system: f64,
}

fn main() -> Result<()> {
    color_eyre::install()?;
    let cli = Cli::parse();

    match cli.command {
        BenchCommand::Check => check_tools(),
        BenchCommand::Fetch(options) => fetch_ena(&options),
        BenchCommand::Prepare(options) => prepare_subsets(&options),
        BenchCommand::Run(options) => run_matrix(&options),
        BenchCommand::Summarize(options) => summarize(&options),
    }
}

fn check_tools() -> Result<()> {
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

fn fetch_ena(options: &CommonOptions) -> Result<()> {
    require_command("curl")?;
    let datasets = read_datasets(&options.config)?;
    let raw_root = options.data_root.join("raw");
    fs::create_dir_all(&raw_root)?;

    for dataset in datasets {
        fetch_dataset(&dataset, &raw_root)?;
    }

    Ok(())
}

fn prepare_subsets(options: &PrepareOptions) -> Result<()> {
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

fn run_matrix(options: &RunOptions) -> Result<()> {
    let paths = ToolPaths::from_environment()?;
    let datasets = read_datasets(&options.common.config)?;
    let run_id = utc_run_id()?;
    let run_dir = options.runs_root.join(&run_id);
    fs::create_dir_all(run_dir.join("metadata"))?;
    write_run_metadata(&run_dir, &run_id, options)?;
    write_tool_versions(&run_dir, &paths)?;

    for dataset in datasets {
        let subset =
            read_subset_metadata(&options.common.data_root, &dataset.name, options.read_pairs)?;
        for tool in &options.tools {
            run_tool(&paths, options, &run_id, &run_dir, &subset, *tool)?;
        }
    }

    eprintln!("Run artifacts: {}", run_dir.display());
    Ok(())
}

fn summarize(options: &SummarizeOptions) -> Result<()> {
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

    for path in result_files {
        let file = File::open(&path)?;
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
    println!("{}", summary_path.display());
    Ok(())
}

impl ToolPaths {
    fn from_environment() -> Result<Self> {
        let file_env = read_env_files()?;
        Ok(Self {
            pairasm: resolve_pairasm(&file_env)?,
            fastp: resolve_tool("FASTP_BIN", "fastp", &file_env)?,
            bbmerge: resolve_tool_any("BBMERGE_BIN", &["bbmerge", "bbmerge.sh"], &file_env)?,
            vsearch: resolve_tool("VSEARCH_BIN", "vsearch", &file_env)?,
            hyperfine: resolve_tool("HYPERFINE_BIN", "hyperfine", &file_env)?,
        })
    }

    fn path_for(&self, tool: Tool) -> &Path {
        match tool {
            Tool::Pairasm => &self.pairasm,
            Tool::Fastp => &self.fastp,
            Tool::Bbmerge => &self.bbmerge,
            Tool::Vsearch => &self.vsearch,
        }
    }
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

fn require_command(binary: &str) -> Result<()> {
    find_on_path(OsStr::new(binary))
        .map(|_| ())
        .ok_or_else(|| color_eyre::eyre::eyre!("required command not found on PATH: {binary}"))
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

fn read_datasets(path: &Path) -> Result<Vec<Dataset>> {
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

fn fetch_dataset(dataset: &Dataset, raw_root: &Path) -> Result<()> {
    let out_dir = raw_root.join(&dataset.name);
    fs::create_dir_all(&out_dir)?;
    let report_path = out_dir.join("ena_filereport.tsv");
    let report_url = format!(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=run_accession,instrument_platform,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_md5&format=tsv&download=true",
        dataset.accession
    );

    eprintln!(
        "Resolving ENA FASTQs for {} ({})",
        dataset.name, dataset.accession
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
        eprintln!("Using cached {}", path.display());
        return Ok(());
    }
    eprintln!("Downloading {url}");
    run_command(Command::new("curl").args(["-fL", url, "-o"]).arg(path))
}

fn validate_gzip(path: &Path) -> Result<()> {
    let mut decoder = GzDecoder::new(File::open(path)?);
    io::copy(&mut decoder, &mut io::sink())?;
    Ok(())
}

fn first_tsv_data_row(path: &Path) -> Result<Vec<String>> {
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

fn read_source_metadata(data_root: &Path, name: &str) -> Result<SourceMetadata> {
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

fn prepare_subset(source: &SourceMetadata, read_pairs: usize, data_root: &Path) -> Result<()> {
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

fn write_first_fastq_records(src: &Path, dst: &Path, read_pairs: usize) -> Result<()> {
    let input = File::open(src)?;
    let decoder = GzDecoder::new(input);
    let mut reader = BufReader::new(decoder);
    let output = File::create(dst)?;
    let encoder = GzEncoder::new(output, Compression::default());
    let mut writer = BufWriter::new(encoder);
    let mut line = String::new();
    let mut lines_written = 0usize;
    let max_lines = read_pairs.saturating_mul(4);

    while lines_written < max_lines {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        writer.write_all(line.as_bytes())?;
        lines_written += 1;
    }
    writer.flush()?;
    Ok(())
}

fn read_subset_metadata(data_root: &Path, name: &str, read_pairs: usize) -> Result<SubsetMetadata> {
    let path = data_root
        .join("subset")
        .join(name)
        .join(format!("{read_pairs}_pairs"))
        .join("subset.tsv");
    let row = first_tsv_data_row(&path)?;
    if row.len() < 5 {
        bail!(
            "subset metadata row had too few columns: {}",
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

fn run_tool(
    paths: &ToolPaths,
    options: &RunOptions,
    run_id: &str,
    run_dir: &Path,
    subset: &SubsetMetadata,
    tool: Tool,
) -> Result<()> {
    let out_dir = run_dir
        .join(&subset.name)
        .join(format!("{}_pairs", subset.read_pairs))
        .join(tool.name());
    fs::create_dir_all(&out_dir)?;
    let command = build_tool_command(paths, options, subset, tool, &out_dir);
    let stdout_log = out_dir.join(format!("{}.stdout.log", tool.name()));
    let stderr_log = out_dir.join(format!("{}.stderr.log", tool.name()));
    let command_string = format!(
        "{} > {} 2> {}",
        shell_join(&command.args),
        shell_quote(&stdout_log.to_string_lossy()),
        shell_quote(&stderr_log.to_string_lossy())
    );
    fs::write(out_dir.join("command.sh"), format!("{command_string}\n"))?;

    eprintln!("[{run_id}] {} {}", subset.name, tool.name());
    let hyperfine_json = out_dir.join("hyperfine.json");
    let hyperfine_md = out_dir.join("hyperfine.md");
    run_command(
        Command::new(&paths.hyperfine)
            .arg("--runs")
            .arg(options.replicates.to_string())
            .arg("--warmup")
            .arg("1")
            .arg("--export-json")
            .arg(&hyperfine_json)
            .arg("--export-markdown")
            .arg(&hyperfine_md)
            .arg("--command-name")
            .arg(tool.name())
            .arg(&command_string),
    )?;

    let report: HyperfineReport = serde_json::from_reader(File::open(&hyperfine_json)?)?;
    let result = report
        .results
        .first()
        .ok_or_else(|| color_eyre::eyre::eyre!("hyperfine report had no results"))?;
    let merged_reads = command
        .merged_output
        .exists()
        .then(|| fastq_record_count(&command.merged_output))
        .transpose()?
        .unwrap_or(0);

    let mut writer = BufWriter::new(File::create(out_dir.join("result.tsv"))?);
    writeln!(
        writer,
        "run_id\tdataset\taccession\tread_pairs\ttool\treplicates\tthreads\toutput_compression\tmean_s\tmedian_s\tstddev_s\tmin_s\tmax_s\tuser_s\tsystem_s\tmerged_reads\tr1_bytes\tr2_bytes\toutput_dir"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        run_id,
        subset.name,
        subset.accession,
        subset.read_pairs,
        command.tool.name(),
        options.replicates,
        options.threads,
        options.output_compression,
        result.mean,
        result.median,
        optional_f64(result.stddev),
        result.min,
        result.max,
        result.user,
        result.system,
        merged_reads,
        file_size(&subset.r1)?,
        file_size(&subset.r2)?,
        out_dir.display()
    )?;
    writer.flush()?;
    Ok(())
}

fn optional_f64(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn build_tool_command(
    paths: &ToolPaths,
    options: &RunOptions,
    subset: &SubsetMetadata,
    tool: Tool,
    out_dir: &Path,
) -> ToolCommand {
    let merged_output = out_dir.join(merged_output_name(tool, options.output_compression));
    let r1 = subset.r1.to_string_lossy().into_owned();
    let r2 = subset.r2.to_string_lossy().into_owned();
    let merged = merged_output.to_string_lossy().into_owned();
    let binary = paths.path_for(tool).to_string_lossy().into_owned();

    let args = match tool {
        Tool::Pairasm => vec![
            binary,
            "-1".to_owned(),
            r1,
            "-2".to_owned(),
            r2,
            "-o".to_owned(),
            merged,
            "--unmerged-out".to_owned(),
            out_dir
                .join("pairasm.unmerged.fastq")
                .to_string_lossy()
                .into_owned(),
            "--summary".to_owned(),
            out_dir
                .join("pairasm.summary.json")
                .to_string_lossy()
                .into_owned(),
            "--progress-every".to_owned(),
            "0".to_owned(),
            "-qqq".to_owned(),
        ],
        Tool::Fastp => vec![
            binary,
            "-i".to_owned(),
            r1,
            "-I".to_owned(),
            r2,
            "--merge".to_owned(),
            "--merged_out".to_owned(),
            merged,
            "--unpaired1".to_owned(),
            out_dir
                .join("fastp.unpaired1.fastq")
                .to_string_lossy()
                .into_owned(),
            "--unpaired2".to_owned(),
            out_dir
                .join("fastp.unpaired2.fastq")
                .to_string_lossy()
                .into_owned(),
            "--failed_out".to_owned(),
            out_dir
                .join("fastp.failed.fastq")
                .to_string_lossy()
                .into_owned(),
            "--thread".to_owned(),
            options.threads.to_string(),
            "--html".to_owned(),
            out_dir.join("fastp.html").to_string_lossy().into_owned(),
            "--json".to_owned(),
            out_dir.join("fastp.json").to_string_lossy().into_owned(),
        ],
        Tool::Bbmerge => vec![
            binary,
            "-da".to_owned(),
            format!("in1={r1}"),
            format!("in2={r2}"),
            format!("out={merged}"),
            format!(
                "outu1={}",
                out_dir.join("bbmerge.unmerged1.fastq").display()
            ),
            format!(
                "outu2={}",
                out_dir.join("bbmerge.unmerged2.fastq").display()
            ),
            format!("threads={}", options.threads),
            "ow=t".to_owned(),
        ],
        Tool::Vsearch => vec![
            binary,
            "--fastq_mergepairs".to_owned(),
            r1,
            "--reverse".to_owned(),
            r2,
            "--fastqout".to_owned(),
            merged,
            "--fastqout_notmerged_fwd".to_owned(),
            out_dir
                .join("vsearch.unmerged1.fastq")
                .to_string_lossy()
                .into_owned(),
            "--fastqout_notmerged_rev".to_owned(),
            out_dir
                .join("vsearch.unmerged2.fastq")
                .to_string_lossy()
                .into_owned(),
            "--threads".to_owned(),
            options.threads.to_string(),
        ],
    };

    ToolCommand {
        tool,
        args,
        merged_output,
    }
}

fn merged_output_name(tool: Tool, output_compression: OutputCompression) -> String {
    match output_compression {
        OutputCompression::Plain => format!("{}.merged.fastq", tool.name()),
        OutputCompression::Gzip => format!("{}.merged.fastq.gz", tool.name()),
    }
}

fn shell_join(args: &[String]) -> String {
    args.iter()
        .map(|arg| shell_quote(arg))
        .collect::<Vec<_>>()
        .join(" ")
}

fn shell_quote(arg: &str) -> String {
    if arg.bytes().all(|byte| {
        byte.is_ascii_alphanumeric()
            || matches!(byte, b'/' | b'.' | b'_' | b'-' | b'=' | b':' | b',')
    }) {
        return arg.to_owned();
    }
    format!("'{}'", arg.replace('\'', "'\\''"))
}

fn write_run_metadata(run_dir: &Path, run_id: &str, options: &RunOptions) -> Result<()> {
    let mut writer = BufWriter::new(File::create(run_dir.join("metadata").join("run.tsv"))?);
    writeln!(writer, "key\tvalue")?;
    writeln!(writer, "run_id\t{run_id}")?;
    writeln!(writer, "read_pairs\t{}", options.read_pairs)?;
    writeln!(writer, "replicates\t{}", options.replicates)?;
    writeln!(writer, "threads\t{}", options.threads)?;
    writeln!(writer, "output_compression\t{}", options.output_compression)?;
    writeln!(writer, "config\t{}", options.common.config.display())?;
    writer.flush()?;
    Ok(())
}

fn write_tool_versions(run_dir: &Path, paths: &ToolPaths) -> Result<()> {
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

fn version_string(path: &Path) -> Result<String> {
    let output = Command::new(path).arg("--version").output()?;
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_owned();
    if !stdout.is_empty() {
        return Ok(stdout);
    }
    Ok(String::from_utf8_lossy(&output.stderr).trim().to_owned())
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

fn fastq_record_count(path: &Path) -> Result<usize> {
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|extension| extension == "gz") {
        Box::new(BufReader::new(GzDecoder::new(File::open(path)?)))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };

    let mut lines = 0usize;
    for line in reader.lines() {
        let _ = line?;
        lines += 1;
    }
    Ok(lines / 4)
}

fn file_size(path: &Path) -> Result<u64> {
    Ok(path.metadata()?.len())
}

fn run_command(command: &mut Command) -> Result<()> {
    let status = command.status()?;
    if !status.success() {
        bail!("command failed with {status}: {:?}", command);
    }
    Ok(())
}

fn utc_run_id() -> Result<String> {
    let output = Command::new("date")
        .args(["-u", "+%Y%m%dT%H%M%SZ"])
        .output()?;
    if !output.status.success() {
        bail!("failed to create UTC run id with date command");
    }
    Ok(String::from_utf8(output.stdout)?.trim().to_owned())
}
