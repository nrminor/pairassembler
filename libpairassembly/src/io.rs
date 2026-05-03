// ROADMAP:
// For I/O to be useful in this library, I'll need to support two or three FASTQ I/O crates, each
// of which must implement some means of extracting the read ID. To make this work, I'll need to
// do the following:
//
// 1. Use noodles by default, making that the first feature, and making noodles installed with the crate by default.
// 2. Make a trait, e.g., `PairAssemble` or `PairOverlap` or something like that, that exposes an `.id()` method (and maybe some other stuff? Maybe an extension trait?) for custom `Record` types or `NewType` pattern wrappers of pre-existing libraries' record types.
// 3. Write a derive macro for the above trait, which can only be derived on structs that have an `id()` method, or something like that.
// 4. Write a bunch of `From<>` and `AsRef<>` impls to make `libpairassembly` useable with noodles and whatever else I decide to support.
// 5. Decide whether I should leave the actual reader and writer stuff to the importing crate.
// 6. Pairing up read mates could be as simple as a blanket implementation of `PartialEq` on all types that implement `PairOverlap`. This could also be a derive macro that is usable when a type first implements `PairOverlap`.

// TODO:
// potential trait to implement to open up the type system:
//
//
// /// Trait for adapting external FASTQ record types into `SequenceRead`.
// pub trait RecordAdapter<'a> {
//     /// Returns the read ID.
//     fn id(&'a self) -> &'a str;

//     /// Returns the sequence as a UTF-8 string.
//     fn sequence(&'a self) -> &'a str;

//     /// Returns the quality scores as a UTF-8 string.
//     fn quality_scores(&'a self) -> &'a str;

//     /// Converts into a `SequenceRead` used internally by `libpairassembly`.
//     fn to_sequence_read(&'a self) -> Result<SequenceRead<'a>> {
//         SequenceRead::try_new(self.id(), self.sequence(), self.quality_scores())
//     }
// }

#[cfg(feature = "noodles")]
pub use noodles::{RecordAdapter, merge_pairs};

#[cfg(feature = "noodles")]
pub mod noodles {
    use crate::{
        Assembler, OverlapParams, OverlapValidator, PairInput, Result, read::SequenceRead,
    };
    use fastq::{
        Record,
        record::{self, Definition},
    };
    use futures::{Stream, StreamExt};
    use noodles_fastq as fastq;
    use std::{borrow::Cow, result::Result as StdResult, str};

    /// Adapter for converting `noodles_fastq::Record` into `SequenceRead`
    #[derive(Debug)]
    pub struct RecordAdapter<'a> {
        pub record: &'a Record,
    }

    impl<'a> From<&'a Record> for SequenceRead<'a> {
        fn from(record: &'a Record) -> Self {
            let id = str::from_utf8(record.name()).unwrap_or("");
            let seq = str::from_utf8(record.sequence()).unwrap_or("");
            let qual = str::from_utf8(record.quality_scores()).unwrap_or("");
            SequenceRead::new(id, seq, qual)
        }
    }

    impl<'a> From<&'a Record> for RecordAdapter<'a> {
        fn from(record: &'a Record) -> Self {
            RecordAdapter { record }
        }
    }

    impl<'a> From<SequenceRead<'a>> for Record {
        fn from(read: SequenceRead<'a>) -> Self {
            let definition = Definition::new(read.id(), "");
            Record::new(definition, read.sequence(), read.quality_scores())
        }
    }

    impl AsRef<Record> for RecordAdapter<'_> {
        fn as_ref(&self) -> &Record {
            self.record
        }
    }

    impl<'a> From<RecordAdapter<'a>> for Record {
        fn from(adapter: RecordAdapter<'a>) -> Self {
            Record::new(
                Definition::new(adapter.record.name(), adapter.record.description()),
                adapter.record.sequence(),
                adapter.record.quality_scores(),
            )
        }
    }

    /// Merges an iterator of paired [`noodles_fastq::Record`]s into validated and corrected records.
    ///
    /// This function takes an iterator over `Result<(Record, Record)>`—typically representing
    /// paired-end reads with matching IDs—and merges them into high-confidence consensus reads using
    /// overlap detection, validation, and base-call correction logic defined in `libpairassembly`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fastq as fastq;
    /// use fastq::{Reader, Record};
    /// use libpairassembly::io::noodles::merge_stream;
    ///
    /// use std::{fs::File, io::BufReader};
    ///
    /// // Open a file containing interleaved paired-end reads
    /// let file = File::open("reads.fastq").expect("failed to open file");
    /// let reader = Reader::new(BufReader::new(file));
    ///
    /// // Read and pair up interleaved FASTQ records
    /// let mut records = reader.records();
    /// let read_pairs = std::iter::from_fn(|| {
    ///     let r1 = records.next()?;
    ///     let r2 = records.next()?;
    ///     Some(r1.and_then(|a| r2.map(|b| (a, b))))
    /// });
    ///
    /// // Merge the pairs into consensus reads
    /// for result in merge_stream(read_pairs) {
    ///     match result {
    ///         Ok(Some(merged)) => {
    ///             println!(">{}", String::from_utf8_lossy(merged.name()));
    ///             println!("{}", String::from_utf8_lossy(merged.sequence()));
    ///             println!("+");
    ///             println!("{}", String::from_utf8_lossy(merged.quality_scores()));
    ///         }
    ///         Ok(None) => {}
    ///         Err(e) => eprintln!("Merge error: {e}"),
    ///     }
    /// }
    /// ```
    ///
    /// # Notes
    ///
    /// - Each pair must have matching IDs to be merged. If IDs do not match, an error will be returned.
    /// - Records are validated and scored according to overlap complexity, expected errors, and mismatch rates.
    /// - Internally, the `merge_stream` function wraps your `Record`s in a lightweight `SequenceRead`
    ///   abstraction and performs fluent method chaining:
    ///
    /// ```rust,ignore
    /// ReadPair::from(read1, read2)?
    ///     .try_find_overlap(params)?
    ///     .try_validate(validator)?
    ///     .merge()
    ///     .correct_quality_scores()
    /// ```
    ///
    /// - This modular pipeline allows advanced users to customize each step of the merging process if needed.
    ///
    pub fn merge_pairs<'a, I>(reads: I) -> impl Iterator<Item = Result<Option<fastq::Record>>> + 'a
    where
        I: Iterator<Item = Result<(Record, Record)>> + 'a,
    {
        reads.filter_map(StdResult::ok).map(handle_pair)
    }

    fn handle_pair(pair: (Record, Record)) -> Result<Option<Record>> {
        let (fwd, rev) = pair;

        let read1 = SequenceRead::from(&fwd);
        let read2 = SequenceRead::from(&rev);

        let pair_input = PairInput::new(read1, read2);

        // Initialize settings for overlapping and for validating those overlaps. We'll just use
        // defaults for demonstration purposes. Note that these are currently consumed, though this
        // may change in the future.
        let overlap_settings = OverlapParams::default();
        let validator = OverlapValidator::default();

        // Use the top-level checked assembler path.
        let assembler = Assembler::builder()
            .overlap(overlap_settings)
            .validate(validator)
            .build()?;

        let Some(merged) = assembler.process_pair(&pair_input)? else {
            return Ok(None);
        };

        let defline = record::Definition::new(merged.id(), "");
        let final_record = Record::new(defline, merged.sequence_bytes(), merged.quality_bytes());

        Ok(Some(final_record))
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "noodles")]
    use crate::{
        io::noodles::{RecordAdapter, merge_pairs},
        prelude::utils::reverse_complement,
    };
    #[cfg(feature = "noodles")]
    use noodles_fastq::record::{Definition, Record};

    #[cfg(feature = "noodles")]
    fn dummy_record(id: &str, seq: &str, qual: &str) -> Record {
        Record::new(Definition::new(id, ""), seq, qual)
    }

    #[test]
    #[cfg(feature = "noodles")]
    fn test_record_adapter_conversion_roundtrip() {
        let record = dummy_record("read1", "ACGT", "IIII");
        let adapter = RecordAdapter::from(&record);
        let roundtrip = Record::from(adapter);

        assert_eq!(record.name(), roundtrip.name());
        assert_eq!(record.sequence(), roundtrip.sequence());
        assert_eq!(record.quality_scores(), roundtrip.quality_scores());
    }

    #[test]
    #[cfg(feature = "noodles")]
    fn test_merge_stream_with_perfect_overlap() {
        let seq =
            "ACGTTGCAGATCTGACCTGAATCGTACGAGTCTAGCGTATGCTAGTCGATCGTACCTGATCGAATCGTAGCTAGTACGATCG";
        let qual = "I".repeat(seq.len());
        let fwd = dummy_record("readX", seq, &qual);
        let rev_seq = reverse_complement(seq);
        let rev = dummy_record("readX", &rev_seq, &qual);

        let pairs = vec![Ok((fwd, rev))];
        let results = merge_pairs(pairs.into_iter()).collect::<Vec<_>>();

        assert_eq!(results.len(), 1);
        let Ok(Some(merged_record)) = &results[0] else {
            panic!("Expected successful merge, got error: {:?}", &results[0]);
        };

        assert_eq!(merged_record.name(), &b"readX"[..]);
        assert_eq!(merged_record.sequence(), seq.as_bytes());
        assert_eq!(merged_record.quality_scores().len(), seq.len());
        assert!(merged_record.quality_scores().iter().all(|q| *q == b'I'));
    }

    #[test]
    #[cfg(feature = "noodles")]
    fn test_merge_stream_mismatched_ids_yields_error() {
        let fwd = dummy_record("read1", "ACGT", "IIII");
        let rev = dummy_record("read2", "TGCA", "IIII");

        let pairs = vec![Ok((fwd, rev))];
        let results = merge_pairs(pairs.into_iter()).collect::<Vec<_>>();

        assert_eq!(results.len(), 1);
        assert!(
            results[0].is_err(),
            "Expected an error due to mismatched IDs"
        );
    }
}
