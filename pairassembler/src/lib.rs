// dev allowances
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]
//
// crate-level lints
#![warn(
    // clippy::pedantic,
    clippy::perf,
    // clippy::todo,
    clippy::unwrap_used,
    clippy::complexity,
    clippy::correctness,
    clippy::absolute_paths,
    clippy::style
)]

use async_compression::tokio::bufread::GzipDecoder;
use color_eyre::eyre::Result;
use futures::TryStreamExt;
use libpairassembly::{BaseCallValidator, OverlapParams, io::merge_pairs};
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
    validation_settings: BaseCallValidator,
    // pub correction_settings:
}

impl RunSettings {
    #[must_use]
    pub fn new(
        overlap_diff_max: usize,
        min_overlap: usize,
        diff_percent_max: f32,
        min_comparisons: usize,
        k: usize,
        min_entropy: usize,
        no_correct: bool,
    ) -> Self {
        // build the overlap settings
        let overlap_settings = OverlapParams::default()
            .with_overlap_diff_max(overlap_diff_max)
            .with_min_overlap(min_overlap)
            .with_diff_percent_max(diff_percent_max)
            .with_min_comparisons(min_comparisons);

        // build the validation settings
        let validation_settings = BaseCallValidator::default()
            .with_k(k)
            .with_min_entropy(min_entropy);

        // correction settings will eventually be built here

        // return the new settings bundle
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
    // Get the path, open it, buffer bytes from the file, and check for a magic number at the top
    // to see if it's gzipped.
    let path = path.as_ref();
    let file_handle = File::open(path).await?;
    let mut read_buffer = BufReader::new(file_handle);
    let is_gzipped = is_gzip_file(&mut read_buffer).await?;

    // pull in the final reader as an owned trait object
    let reader: Box<dyn AsyncBufRead + Unpin + Send> =
        // box new buffered readers on the heap based on whether the file needs to be decoded or not
        if is_gzipped {
            Box::new(BufReader::new(GzipDecoder::new(read_buffer)))
        } else {
            Box::new(read_buffer)
        };

    Ok(AsyncReader::new(reader))
}

#[allow(clippy::absolute_paths)]
async fn open_fastq_reader(
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
    use std::{any::Any, future, path::PathBuf, sync::Arc};

    use futures::StreamExt;
    use libpairassembly::{ReadPair, SequenceRead, errors::PairingError::UnmatchedIds};
    use noodles::fastq::Record as FastqRecord;
    use tokio::task;

    use super::*;

    /// Enum for storing information about whether a read could be merged or not
    enum OverlapResult<T> {
        Overlap(T),
        NoOverlap((T, T)),
    }
    use OverlapResult::*;

    pub async fn run(
        input1: String,
        input2: Option<String>,
        output_file: Option<String>,
        unmerged_output: Option<String>,
        settings: RunSettings,
    ) -> color_eyre::Result<()> {
        let mut fastq_reader1 = open_async_fastq_reader(&input1).await?;
        let Some(input2) = input2 else {
            // TODO: Need to write an implementation that converts interleaved pairs into an iterator of mate tuples
            unimplemented!()
        };
        let mut fastq_reader2 = open_async_fastq_reader(&input2).await?;

        let merged_reads = fastq_reader1
            .records()
            .zip(fastq_reader2.records())
            .filter_map(|(attempt1, attempt2)| async move {
                match (attempt1, attempt2) {
                    (Ok(fwd), Ok(rev)) => Some((fwd, rev)),
                    _ => None, // TODO: Errors should not be ignored obvi
                }
            })
            .scan(0_u8, |id_mismatch_tally, pair| {
                // get references to each record in the pair and increment the mismatch tally if the
                // names don't match as is expected in paired FASTQs
                let (ref fwd, ref rev) = pair;
                if fwd.name() != rev.name() {
                    *id_mismatch_tally += 1;
                }

                // put the resulting pair in a ready future, which, because it does not need to be
                // awaited and polled, will spare us most of the complexity of managing lifetimes
                // in async contexts. NOTE: This is a really nice trick and could be helpful in
                // a variety of contexts where tokio is used to do compute-bound as well as IO-bound
                // work.
                future::ready(if *id_mismatch_tally < 3 {
                    Some(pair)
                } else {
                    // TODO: We probably need some error handling here
                    None
                })
            })
            // By this point in the method chain, we have zipped together read mates and removed records
            // that produced errors. We've also scanned for up to three mismatching read IDs, in which
            // case the the iterator prematurely terminates. If all is well, we can proceed by taking
            // ownership of the read data, bundle them into a `ReadMates` instance, and run the
            // pairassembly workflow on it.
            .then(move |(fwd, rev)| {
                task::spawn_blocking(
                    move || -> libpairassembly::Result<OverlapResult<FastqRecord>> {
                        // read id handling. This is no longer necessary here and is in fact redundant.
                        let fwd_id = fwd.definition().name();
                        let rev_id = rev.definition().name();
                        if fwd_id != rev_id {
                            return Err(UnmatchedIds(fwd_id.to_string(), rev_id.to_string()).into());
                        }

                        // repackage reads into mates. Note that this does not involve any copying of
                        // data and is instead merely to help the compiler ensure what we're doing is
                        // correct.
                        let read1 = SequenceRead::from(&fwd);
                        let read2 = SequenceRead::from(&rev);
                        let mates = ReadPair::from(read1, read2)?;

                        // Initialize settings for overlapping and for validating those overlaps. We'll just use
                        // defaults for demonstration purposes. Note that these are currently consumed, though this
                        // may change in the future. We don't need to clone here because the settings are `Copy`.
                        let overlap_settings = settings.overlap_settings;
                        let validator = settings.validation_settings;

                        // First, search for an overlap, early-returning a `NoOverlap` if there is none
                        let Some(overlap) = mates.overlap(&overlap_settings)? else {
                            let res = NoOverlap((fwd, rev));
                            return Ok(res);
                        };

                        // Same thing with validation -- if there isn't a valid overlap, early-return the
                        // original record
                        let Ok(validated) = overlap.validate(&mates, &validator) else {
                            let res = NoOverlap((fwd, rev));
                            return Ok(res);
                        };

                        // if there is an overlap, proceed to validation and merging. Unlike the above early-
                        // return cases, we should just propogate an error if one occurs here.
                        let merged = validated.merge()?;

                        // skip correction and early-return if not turned off
                        if settings.no_correct {
                            let final_record = FastqRecord::new(
                                fwd.definition().clone(),
                                merged.sequence(),
                                merged.qualities(),
                            );
                            return Ok(Overlap(final_record));
                        }

                        // otherwise, run correction
                        let corrected = merged.correct()?;

                        // make a new record and return
                        let final_record = FastqRecord::new(
                            fwd.definition().clone(),
                            corrected.sequence_bytes(),
                            corrected.quality_bytes(),
                        );

                        Ok(Overlap(final_record))
                    },
                )
            })
            .then(|join_result| async move {
                join_result.unwrap() // join handle unwrap
            });

        todo!()
    }

    #[allow(unused_variables)]
    pub async fn run_sync(
        input1: String,
        input2: Option<String>,
        output_file: Option<String>,
        unmerged_output: Option<String>,
        settings: RunSettings,
    ) -> color_eyre::Result<()> {
        let mut fastq_reader1 = open_fastq_reader(&input1).await?;
        let Some(input2) = input2 else {
            unimplemented!()
        };
        let mut fastq_reader2 = open_fastq_reader(&input2).await?;

        let putative_mates = fastq_reader1
            .records()
            .flat_map(|record| record.ok())
            .zip(fastq_reader2.records().flat_map(|record| record.ok()));

        let mut mismatch_counter = 0_usize;
        for (fwd, rev) in putative_mates {
            if fwd.definition().name() != rev.definition().name() {
                mismatch_counter += 1;
                if mismatch_counter > 3 {
                    break;
                }
            }

            let read1 = SequenceRead::from(&fwd);
            let read2 = SequenceRead::from(&rev);

            let mates = ReadPair::from(read1, read2)?;

            // Initialize settings for overlapping and for validating those overlaps. We'll just use
            // defaults for demonstration purposes. Note that these are currently consumed, though this
            // may change in the future.
            let overlap_settings = settings.overlap_settings;
            let validator = settings.validation_settings;

            // Use a chain of methods to execute the whole pipeline
            let merged = mates
                .overlap(&overlap_settings)?
                .ok_or_else(|| color_eyre::eyre::eyre!("No overlap found."))?
                .validate(&mates, &validator)?
                .merge()?
                .correct()?;

            let final_record = FastqRecord::new(
                fwd.definition().clone(),
                merged.sequence_bytes(),
                merged.quality_bytes(),
            );
        }

        todo!()
    }
}
