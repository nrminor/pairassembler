use std::io::{self, IsTerminal};

use console::style;

use crate::model::Tool;

pub fn heading_stdout(title: &str) -> String {
    if stdout_styled() {
        style(title).cyan().bold().to_string()
    } else {
        title.to_owned()
    }
}

pub fn tool_name_stdout(tool_name: &str) -> String {
    styled_tool_name(tool_name, stdout_styled())
}

pub fn tool_stderr(tool: Tool) -> String {
    styled_tool_name(tool.name(), stderr_styled())
}

pub fn dataset_stderr(dataset_name: &str) -> String {
    if stderr_styled() {
        style(dataset_name).bold().to_string()
    } else {
        dataset_name.to_owned()
    }
}

pub fn muted_stderr(text: impl std::fmt::Display) -> String {
    let text = text.to_string();
    if stderr_styled() {
        style(text).dim().to_string()
    } else {
        text
    }
}

pub fn path_stderr(path: impl std::fmt::Display) -> String {
    let path = path.to_string();
    if stderr_styled() {
        style(path).underlined().to_string()
    } else {
        path
    }
}

pub fn print_workflow_phase(step: &str, title: &str) {
    if stdout_styled() {
        println!();
        println!(
            "{} {}",
            style("━━ Benchmark workflow").cyan().bold(),
            style(step).cyan().bold()
        );
        println!("{}", style(title).dim());
    } else {
        println!();
        println!("━━ Benchmark workflow {step}");
        println!("{title}");
    }
}

pub fn color_tool_names_for_stdout(text: &str) -> String {
    if !stdout_styled() {
        return text.to_owned();
    }

    ["pairasm", "fastp", "bbmerge", "vsearch"]
        .into_iter()
        .fold(text.to_owned(), |text, tool| {
            text.replace(tool, &tool_name_stdout(tool))
        })
}

fn styled_tool_name(tool_name: &str, styled: bool) -> String {
    if !styled {
        return tool_name.to_owned();
    }

    match tool_name {
        "pairasm" => style(tool_name).cyan().bold().to_string(),
        "fastp" => style(tool_name).green().bold().to_string(),
        "bbmerge" => style(tool_name).magenta().bold().to_string(),
        "vsearch" => style(tool_name).yellow().bold().to_string(),
        _ => style(tool_name).bold().to_string(),
    }
}

fn stdout_styled() -> bool {
    io::stdout().is_terminal()
}

fn stderr_styled() -> bool {
    io::stderr().is_terminal()
}
