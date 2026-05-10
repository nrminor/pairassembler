use clap::ValueEnum;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum Tool {
    Pairasm,
    Fastp,
    Bbmerge,
    Vsearch,
}

impl Tool {
    fn as_str(self) -> &'static str {
        match self {
            Tool::Pairasm => "pairasm",
            Tool::Fastp => "fastp",
            Tool::Bbmerge => "bbmerge",
            Tool::Vsearch => "vsearch",
        }
    }
}

impl std::fmt::Display for Tool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}
