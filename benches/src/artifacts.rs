use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ArtifactKind {
    Command,
    StdoutLog,
    StderrLog,
    HyperfineJson,
    MergedFastq,
    PairasmSummaryJson,
    PairasmUnmergedFastq,
    FastpJson,
    FastpHtml,
    FastpUnpaired1Fastq,
    FastpUnpaired2Fastq,
    FastpFailedFastq,
    BbmergeUnmerged1Fastq,
    BbmergeUnmerged2Fastq,
    VsearchUnmerged1Fastq,
    VsearchUnmerged2Fastq,
}

impl ArtifactKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Command => "command",
            Self::StdoutLog => "stdout_log",
            Self::StderrLog => "stderr_log",
            Self::HyperfineJson => "hyperfine_json",
            Self::MergedFastq => "merged_fastq",
            Self::PairasmSummaryJson => "pairasm_summary_json",
            Self::PairasmUnmergedFastq => "pairasm_unmerged_fastq",
            Self::FastpJson => "fastp_json",
            Self::FastpHtml => "fastp_html",
            Self::FastpUnpaired1Fastq => "fastp_unpaired1_fastq",
            Self::FastpUnpaired2Fastq => "fastp_unpaired2_fastq",
            Self::FastpFailedFastq => "fastp_failed_fastq",
            Self::BbmergeUnmerged1Fastq => "bbmerge_unmerged1_fastq",
            Self::BbmergeUnmerged2Fastq => "bbmerge_unmerged2_fastq",
            Self::VsearchUnmerged1Fastq => "vsearch_unmerged1_fastq",
            Self::VsearchUnmerged2Fastq => "vsearch_unmerged2_fastq",
        }
    }
}

impl std::fmt::Display for ArtifactKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ArtifactRequirement {
    Required,
    Optional,
}

impl ArtifactRequirement {
    pub(crate) fn is_required(self) -> bool {
        matches!(self, Self::Required)
    }
}

pub(crate) struct ArtifactRecord<'a> {
    pub(crate) kind: ArtifactKind,
    pub(crate) path: &'a Path,
    pub(crate) requirement: ArtifactRequirement,
}
