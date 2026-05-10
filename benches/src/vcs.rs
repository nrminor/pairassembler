use std::process::Command;

#[derive(Debug, Default)]
pub(crate) struct VcsMetadata {
    pub(crate) vcs_kind: Option<String>,
    pub(crate) change_id: Option<String>,
    pub(crate) commit_id: Option<String>,
    pub(crate) description: Option<String>,
    pub(crate) working_copy_dirty: Option<bool>,
}

pub(crate) fn collect_vcs_metadata() -> VcsMetadata {
    collect_jj_metadata()
        .or_else(collect_git_metadata)
        .unwrap_or_default()
}

fn collect_jj_metadata() -> Option<VcsMetadata> {
    let status = command_stdout(Command::new("jj").arg("status"))?;
    let log = command_stdout(Command::new("jj").arg("log").args([
        "-r",
        "@",
        "--no-graph",
        "--limit",
        "1",
        "--template",
        "change_id ++ \"\\n\" ++ commit_id ++ \"\\n\" ++ description.first_line() ++ \"\\n\"",
    ]))?;
    Some(VcsMetadata {
        vcs_kind: Some("jj".to_owned()),
        change_id: parse_jj_change_id(&log),
        commit_id: parse_jj_commit_id(&log),
        description: parse_jj_description(&log),
        working_copy_dirty: Some(status.contains("Working copy changes:")),
    })
}

fn collect_git_metadata() -> Option<VcsMetadata> {
    let commit_id = command_stdout(Command::new("git").args(["rev-parse", "HEAD"]))?;
    let description = command_stdout(Command::new("git").args(["log", "-1", "--format=%s"]))?;
    let status = command_stdout(Command::new("git").args(["status", "--porcelain"]))?;
    Some(VcsMetadata {
        vcs_kind: Some("git".to_owned()),
        change_id: None,
        commit_id: Some(commit_id.trim().to_owned()),
        description: Some(description.trim().to_owned()),
        working_copy_dirty: Some(!status.trim().is_empty()),
    })
}

fn command_stdout(command: &mut Command) -> Option<String> {
    let output = command.output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8(output.stdout).ok())?
}

fn parse_jj_change_id(log: &str) -> Option<String> {
    log.lines()
        .next()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
}

fn parse_jj_commit_id(log: &str) -> Option<String> {
    log.lines()
        .nth(1)
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
}

fn parse_jj_description(log: &str) -> Option<String> {
    log.lines()
        .nth(2)
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
}
