use std::io;

use color_eyre::eyre::{Result, bail};
use tabled::{Table, Tabled, settings::Style};

use crate::{
    cli::{ReportCommand, ReportOptions, RunScopedReportOptions},
    db::{AgreementRow, BenchmarkDb, RunSummary, ToolResultRow},
    ui,
};

pub fn report(options: &ReportOptions) -> Result<()> {
    match &options.command {
        ReportCommand::ReadIdOverlap(report_options) => report_read_id_overlap(report_options),
        ReportCommand::ToolResultsTsv(report_options) => report_tool_results_tsv(report_options),
        ReportCommand::TimingMarkdown(report_options) => report_timing_markdown(report_options),
    }
}

fn report_read_id_overlap(options: &RunScopedReportOptions) -> Result<()> {
    let database = BenchmarkDb::open_existing(&options.db)?;
    let run_key = report_run_key(&database, options)?;
    let summary = database.run_summary(&run_key)?;
    let rows = database.agreement_rows(&run_key)?;
    write_read_id_overlap_report(&summary, &rows)
}

fn report_tool_results_tsv(options: &RunScopedReportOptions) -> Result<()> {
    let database = BenchmarkDb::open_existing(&options.db)?;
    let run_key = report_run_key(&database, options)?;
    let rows = database.tool_result_rows(&run_key)?;
    write_tool_results_tsv(&run_key, &rows)
}

fn report_timing_markdown(options: &RunScopedReportOptions) -> Result<()> {
    let database = BenchmarkDb::open_existing(&options.db)?;
    let run_key = report_run_key(&database, options)?;
    let summary = database.run_summary(&run_key)?;
    let rows = database.tool_result_rows(&run_key)?;
    write_timing_markdown(&summary, &rows)
}

fn report_run_key(database: &BenchmarkDb, options: &RunScopedReportOptions) -> Result<String> {
    options
        .run
        .clone()
        .map(Ok)
        .unwrap_or_else(|| database.latest_completed_run_key_for_mode(&options.mode.to_string()))
}

fn write_read_id_overlap_report(summary: &RunSummary, rows: &[AgreementRow]) -> Result<()> {
    if rows.is_empty() {
        bail!(
            "no merged read-ID set overlap report is available for run {}\n\nOverlap reporting needs at least two tool executions in the same benchmark run.\n\nTo run the standard real-data comparison:\n\n  just benchmark\n\nFor tuned/comparability mode:\n\n  just benchmark-tuned",
            summary.run_key
        );
    }

    print_run_header("Merged read-ID set overlap", summary);
    let rows = rows
        .iter()
        .map(AgreementReportRow::from)
        .collect::<Vec<_>>();
    let mut table = Table::new(rows);
    table.with(Style::rounded());
    println!("{}", ui::color_tool_names_for_stdout(&table.to_string()));

    Ok(())
}

fn write_tool_results_tsv(run_label: &str, rows: &[ToolResultRow]) -> Result<()> {
    if rows.is_empty() {
        bail!("no tool results are available for run {run_label}");
    }

    let stdout = io::stdout();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(stdout.lock());

    writer.write_record([
        "run_id",
        "benchmark_mode",
        "dataset",
        "accession",
        "read_pairs",
        "tool",
        "replicates",
        "threads",
        "mean_s",
        "median_s",
        "stddev_s",
        "min_s",
        "max_s",
        "user_s",
        "system_s",
        "merged_reads",
        "merged_read_records",
        "r1_bytes",
        "r2_bytes",
        "output_dir",
    ])?;

    for row in rows {
        writer.write_record([
            row.run_key.clone(),
            row.benchmark_mode.clone(),
            row.dataset_name.clone(),
            row.accession.clone(),
            row.read_pairs.to_string(),
            row.tool.clone(),
            row.replicates.to_string(),
            row.threads.to_string(),
            row.mean_s.to_string(),
            row.median_s.to_string(),
            optional_f64(row.stddev_s),
            row.min_s.to_string(),
            row.max_s.to_string(),
            row.user_s.to_string(),
            row.system_s.to_string(),
            row.merged_reads.to_string(),
            row.merged_read_records.to_string(),
            row.r1_bytes.to_string(),
            row.r2_bytes.to_string(),
            row.output_dir.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_timing_markdown(summary: &RunSummary, rows: &[ToolResultRow]) -> Result<()> {
    if rows.is_empty() {
        bail!("no timing rows are available for run {}", summary.run_key);
    }

    print_run_header("Tool timing summary", summary);
    println!("| dataset | tool | mean_s | median_s | stddev_s | min_s | max_s | merged_reads |");
    println!("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |");
    for row in rows {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} |",
            row.dataset_name,
            row.tool,
            row.mean_s,
            row.median_s,
            optional_f64(row.stddev_s),
            row.min_s,
            row.max_s,
            row.merged_reads,
        );
    }
    Ok(())
}

fn print_run_header(title: &str, summary: &RunSummary) {
    println!("{}", ui::heading_stdout(title));
    println!("run: {}", summary.run_key);
    println!("mode: {}", summary.benchmark_mode);
    println!("read_pairs: {}", summary.read_pairs);
    println!("replicates: {}", summary.replicates);
    println!("threads: {}", summary.threads);
    println!("status: {}", summary.status);
    if let Some(completed_at) = &summary.completed_at {
        println!("completed_at: {completed_at}");
    }
    println!("run_dir: {}", summary.run_dir);
    println!();
}

fn optional_f64(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

#[derive(Tabled)]
struct AgreementReportRow<'a> {
    dataset: &'a str,
    left_tool: &'a str,
    right_tool: &'a str,
    left_merged: i64,
    right_merged: i64,
    shared_merged: i64,
    left_only: i64,
    right_only: i64,
    jaccard: f64,
}

impl<'a> From<&'a AgreementRow> for AgreementReportRow<'a> {
    fn from(row: &'a AgreementRow) -> Self {
        Self {
            dataset: &row.dataset_name,
            left_tool: &row.left_tool,
            right_tool: &row.right_tool,
            left_merged: row.left_merged,
            right_merged: row.right_merged,
            shared_merged: row.shared_merged,
            left_only: row.left_only,
            right_only: row.right_only,
            jaccard: row.jaccard,
        }
    }
}
