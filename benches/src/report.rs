use color_eyre::eyre::{Result, bail};
use tabled::{Table, Tabled, settings::Style};

use crate::{
    cli::{ReportCommand, ReportOptions},
    db::{AgreementRow, BenchmarkDb},
};

pub fn report(options: &ReportOptions) -> Result<()> {
    match &options.command {
        ReportCommand::Agreement(agreement_options) => {
            let database = BenchmarkDb::open(&agreement_options.db)?;
            let run_label = match &agreement_options.run {
                Some(run_label) => run_label.clone(),
                None => database.latest_run_label()?,
            };
            let rows = database.agreement_rows(&run_label)?;
            write_agreement_report(&run_label, &rows)
        },
    }
}

fn write_agreement_report(run_label: &str, rows: &[AgreementRow]) -> Result<()> {
    if rows.is_empty() {
        bail!(
            "no tool-agreement report is available for run {run_label}\n\nAgreement reporting needs merged reads from at least two tools in the same benchmark run.\n\nTo run the standard real-data comparison:\n\n  just benchmark\n\nFor tuned/comparability mode:\n\n  just benchmark-tuned"
        );
    }

    let rows = rows
        .iter()
        .map(AgreementReportRow::from)
        .collect::<Vec<_>>();
    let mut table = Table::new(rows);
    table.with(Style::rounded());
    println!("{table}");

    Ok(())
}

#[derive(Tabled)]
struct AgreementReportRow<'a> {
    #[tabled(rename = "run")]
    run_label: &'a str,
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
            run_label: &row.run_label,
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
