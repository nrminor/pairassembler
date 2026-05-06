use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use color_eyre::eyre::{Result, WrapErr};
use flate2::{Compression, write::GzEncoder};
use libpairassembly::OwnedSequenceRead;
use noodles::fastq::{Record as FastqRecord, io::Writer as FastqWriter, record::Definition};

use crate::RunRequest;

use super::is_gzip_path;

enum FastqOutput {
    Plain(FastqWriter<Box<dyn Write + Send>>),
    Gzip(FastqWriter<GzEncoder<Box<dyn Write + Send>>>),
}

impl FastqOutput {
    fn new(path: Option<&str>) -> Result<Self> {
        match path {
            Some(path) if is_gzip_path(Path::new(path)) => {
                let file = File::create(path).wrap_err_with(|| {
                    format!(
                        "failed to create gzip FASTQ output: {path}\nhelp: check that the parent directory exists and is writable"
                    )
                })?;
                let writer: Box<dyn Write + Send> = Box::new(BufWriter::new(file));
                Ok(Self::Gzip(FastqWriter::new(GzEncoder::new(
                    writer,
                    Compression::default(),
                ))))
            },
            Some(path) => {
                let file = File::create(path).wrap_err_with(|| {
                    format!(
                        "failed to create FASTQ output: {path}\nhelp: check that the parent directory exists and is writable"
                    )
                })?;
                Ok(Self::Plain(FastqWriter::new(Box::new(BufWriter::new(
                    file,
                )))))
            },
            None => Ok(Self::Plain(FastqWriter::new(Box::new(BufWriter::new(
                io::stdout(),
            ))))),
        }
    }

    fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.write_record(record)?,
            Self::Gzip(writer) => writer.write_record(record)?,
        }
        Ok(())
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => writer.get_mut().flush()?,
            Self::Gzip(writer) => {
                let encoder = writer.into_inner();
                let mut inner = encoder.finish()?;
                inner.flush()?;
            },
        }
        Ok(())
    }
}

pub(super) struct OutputHandles {
    merged: FastqOutput,
    unmerged: Option<FastqOutput>,
}

impl OutputHandles {
    pub(super) fn new(request: &RunRequest) -> Result<Self> {
        Ok(Self {
            merged: FastqOutput::new(request.output_file.as_deref())?,
            unmerged: match request.unmerged_output.as_deref() {
                Some(path) => Some(FastqOutput::new(Some(path))?),
                None => None,
            },
        })
    }

    pub(super) fn write_merged(&mut self, read: &OwnedSequenceRead) -> Result<()> {
        self.merged
            .write_record(&merged_record(read))
            .wrap_err_with(|| {
                format!(
                    "failed to write merged FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                    read.id()
                )
            })
    }

    pub(super) fn write_unmerged_records(
        &mut self,
        r1: &FastqRecord,
        r2: &FastqRecord,
    ) -> Result<bool> {
        let Some(output) = &mut self.unmerged else {
            return Ok(false);
        };

        output.write_record(r1).wrap_err_with(|| {
            format!(
                "failed to write unmerged R1 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(r1.name().as_ref())
            )
        })?;
        output.write_record(r2).wrap_err_with(|| {
            format!(
                "failed to write unmerged R2 FASTQ record\nread_id: {}\nhelp: check downstream pipe or output filesystem health",
                String::from_utf8_lossy(r2.name().as_ref())
            )
        })?;

        Ok(true)
    }

    pub(super) fn finish(self) -> Result<()> {
        self.merged
            .finish()
            .wrap_err("failed to finalize merged FASTQ output")?;
        if let Some(output) = self.unmerged {
            output
                .finish()
                .wrap_err("failed to finalize unmerged FASTQ output")?;
        }
        Ok(())
    }
}

fn merged_record(read: &OwnedSequenceRead) -> FastqRecord {
    FastqRecord::new(
        Definition::new(read.id(), ""),
        read.sequence_bytes(),
        read.quality_bytes(),
    )
}
