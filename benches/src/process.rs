use std::process::Command;

use color_eyre::eyre::{Result, bail};

pub fn run_command(command: &mut Command) -> Result<()> {
    let status = command.status()?;
    if !status.success() {
        bail!("command failed with {status}: {:?}", command);
    }
    Ok(())
}
