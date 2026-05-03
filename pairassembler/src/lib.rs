#![allow(dead_code, unused_imports, unused_variables, unused_mut)]
#![warn(
    clippy::perf,
    clippy::unwrap_used,
    clippy::complexity,
    clippy::correctness,
    clippy::absolute_paths,
    clippy::style
)]

use async_compression::tokio::bufread::GzipDecoder;
use color_eyre::eyre::Result;
use futures::TryStreamExt;
use libpairassembly::{Assembler, OverlapParams, OverlapValidator, PairInput};
use noodles::fastq::AsyncReader;
use std::path::Path;
use tokio::{
    fs::File,
    io::{AsyncBufRead, AsyncBufReadExt, AsyncRead, BufReader},
};

pub mod cli;

#[derive(Debug)]
pub struct RunSettings {
    no_correct: bool,
    overlap_settings: OverlapParams,
    validation_settings: OverlapValidator,
}

impl RunSettings {
    #[must_use]
    pub fn new(
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        k: usize,
        min_complexity_score: usize,
        no_correct: bool,
    ) -> Self {
        let overlap_settings = OverlapParams::default()
            .with_overlap_diff_max(overlap_diff_max)
            .with_min_overlap(min_overlap)
            .with_diff_percent_max(diff_percent_max)
            .with_min_comparisons(min_comparisons);

        let validation_settings = OverlapValidator::default()
            .with_k(k)
            .with_min_complexity_score(min_complexity_score);

        RunSettings {
            no_correct,
            overlap_settings,
            validation_settings,
        }
    }
}

async fn open_async_fastq_reader(
    path: impl AsRef<Path>,
) -> Result<AsyncReader<Box<dyn AsyncBufRead + Unpin + Send>>> {
    let path = path.as_ref();
    let file_handle = File::open(path).await?;
    let mut read_buffer = BufReader::new(file_handle);
    let is_gzipped = is_gzip_file(&mut read_buffer).await?;

    let reader: Box<dyn AsyncBufRead + Unpin + Send> = if is_gzipped {
        Box::new(BufReader::new(GzipDecoder::new(read_buffer)))
    } else {
        Box::new(read_buffer)
    };

    Ok(AsyncReader::new(reader))
}

#[allow(clippy::absolute_paths)]
fn open_fastq_reader(
    path: impl AsRef<Path>,
) -> Result<noodles::fastq::Reader<Box<dyn std::io::BufRead>>> {
    let path = path.as_ref();
    let file_handle = std::fs::File::open(path)?;
    let mut read_buffer = std::io::BufReader::new(file_handle);

    let reader: Box<dyn std::io::BufRead> = if path.ends_with("gz") {
        Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(
            read_buffer,
        )))
    } else {
        Box::new(read_buffer)
    };
    let reader = noodles::fastq::Reader::new(reader);
    Ok(reader)
}

async fn is_gzip_file(file: &mut BufReader<File>) -> Result<bool> {
    let buffer = file.fill_buf().await?;
    Ok(buffer.len() >= 2 && buffer[0] == 0x1F && buffer[1] == 0x8B)
}

pub mod merging {
    use std::{any::Any, future, path::PathBuf, result::Result as StdResult, str, sync::Arc};

    use futures::StreamExt;
    use libpairassembly::{OverlapSearch, SequenceRead, errors::PairingError::UnmatchedIds};
    use noodles::fastq::Record as FastqRecord;
    use tokio::task;

    use super::{Assembler, PairInput, RunSettings, open_async_fastq_reader, open_fastq_reader};

    /// Output branch for a processed read pair.
    enum OverlapResult<T> {
        Overlap(T),
        NoOverlap((T, T)),
    }

    fn sequence_read_from_record(record: &FastqRecord) -> color_eyre::Result<SequenceRead<'_>> {
        let id = str::from_utf8(record.name().as_ref())?;
        let seq = str::from_utf8(record.sequence())?;
        let qual = str::from_utf8(record.quality_scores())?;

        Ok(SequenceRead::try_new(id, seq, qual)?)
    }

    /// Run asynchronous pair merging over two FASTQ inputs.
    ///
    /// # Errors
    ///
    /// Returns an error when input files cannot be read/decoded or when pair
    /// processing fails in overlap, validation, merge, or correction stages.
    ///
    /// # Panics
    ///
    /// Panics if a merge worker task panics while executing in `spawn_blocking`.
    pub async fn run(
        input1: String,
        input2: Option<String>,
        output_file: Option<String>,
        unmerged_output: Option<String>,
        settings: RunSettings,
    ) -> color_eyre::Result<()> {
        let mut fastq_reader1 = open_async_fastq_reader(&input1).await?;
        let Some(input2) = input2 else {
            unimplemented!()
        };
        let mut fastq_reader2 = open_async_fastq_reader(&input2).await?;

        let merged_reads =
            fastq_reader1
                .records()
                .zip(fastq_reader2.records())
                .filter_map(|(attempt1, attempt2)| async move {
                    match (attempt1, attempt2) {
                        (Ok(fwd), Ok(rev)) => Some((fwd, rev)),
                        _ => None,
                    }
                })
                .scan(0_u8, |id_mismatch_tally, pair| {
                    let (ref fwd, ref rev) = pair;
                    if fwd.name() != rev.name() {
                        *id_mismatch_tally += 1;
                    }

                    future::ready(if *id_mismatch_tally < 3 {
                        Some(pair)
                    } else {
                        None
                    })
                })
                .then(move |(fwd, rev)| {
                    task::spawn_blocking(
                        move || -> color_eyre::Result<OverlapResult<FastqRecord>> {
                            let fwd_id = fwd.definition().name();
                            let rev_id = rev.definition().name();
                            if fwd_id != rev_id {
                                return Err(
                                    UnmatchedIds(fwd_id.to_string(), rev_id.to_string()).into()
                                );
                            }

                            let read1 = sequence_read_from_record(&fwd)?;
                            let read2 = sequence_read_from_record(&rev)?;
                            let pair_input = PairInput::new(read1, read2);

                            let overlap_settings = settings.overlap_settings;
                            let validator = settings.validation_settings;

                            let assembler = Assembler::builder()
                                .with_overlap_params(overlap_settings)
                                .with_validator(validator)
                                .build()?;

                            let overlap_ctx =
                                match assembler.on_pair(&pair_input)?.find_overlap()? {
                                    OverlapSearch::Found(ctx) => ctx,
                                    OverlapSearch::NoOverlap(_) => {
                                        return Ok(OverlapResult::NoOverlap((fwd, rev)));
                                    },
                                };

                            let Ok(validated_ctx) = overlap_ctx.validate() else {
                                return Ok(OverlapResult::NoOverlap((fwd, rev)));
                            };

                            let merged = validated_ctx.merge()?;

                            if settings.no_correct {
                                let merged = merged.into_owned_read()?;
                                let final_record = FastqRecord::new(
                                    fwd.definition().clone(),
                                    merged.sequence().as_bytes(),
                                    merged.quality_scores().as_bytes(),
                                );
                                return Ok(OverlapResult::Overlap(final_record));
                            }

                            let corrected = merged.correct()?.into_owned_read()?;

                            let final_record = FastqRecord::new(
                                fwd.definition().clone(),
                                corrected.sequence().as_bytes(),
                                corrected.quality_scores().as_bytes(),
                            );

                            Ok(OverlapResult::Overlap(final_record))
                        },
                    )
                })
                .then(|join_result| async move {
                    join_result.expect("merge worker task should not panic")
                });

        todo!()
    }

    #[allow(unused_variables)]
    #[allow(clippy::needless_pass_by_value)]
    /// Run synchronous pair merging over two FASTQ inputs.
    ///
    /// # Errors
    ///
    /// Returns an error when input files cannot be read/decoded or when pair
    /// processing fails in overlap, validation, merge, or correction stages.
    pub fn run_sync(
        input1: String,
        input2: Option<String>,
        output_file: Option<String>,
        unmerged_output: Option<String>,
        settings: RunSettings,
    ) -> color_eyre::Result<()> {
        let mut fastq_reader1 = open_fastq_reader(&input1)?;
        let Some(input2) = input2 else {
            unimplemented!()
        };
        let mut fastq_reader2 = open_fastq_reader(&input2)?;

        let putative_mates = fastq_reader1
            .records()
            .filter_map(StdResult::ok)
            .zip(fastq_reader2.records().filter_map(StdResult::ok));

        let assembler = Assembler::builder()
            .with_overlap_params(settings.overlap_settings)
            .with_validator(settings.validation_settings)
            .build()?;

        let mut mismatch_counter = 0_usize;
        for (fwd, rev) in putative_mates {
            if fwd.definition().name() != rev.definition().name() {
                mismatch_counter += 1;
                if mismatch_counter > 3 {
                    break;
                }
            }
            let read1 = sequence_read_from_record(&fwd)?;
            let read2 = sequence_read_from_record(&rev)?;

            let pair_input = PairInput::new(read1, read2);
            let search = assembler.on_pair(&pair_input)?.find_overlap()?;
            let OverlapSearch::Found(overlap) = search else {
                continue;
            };
            let merged = overlap.validate()?.merge()?.correct()?.into_owned_read()?;

            let final_record = FastqRecord::new(
                fwd.definition().clone(),
                merged.sequence().as_bytes(),
                merged.quality_scores().as_bytes(),
            );
        }

        todo!()
    }
}
