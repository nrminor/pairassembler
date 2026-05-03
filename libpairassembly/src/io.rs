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

    /// Adapter for converting `noodles_fastq::Record` values into library read views.
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

    /// Merge paired [`noodles_fastq::Record`] values with the default assembler pipeline.
    ///
    /// Each input item is a pair of mate records with matching IDs. A pair that has an acceptable
    /// overlap is returned as `Ok(Some(record))`; a successfully processed pair with no acceptable
    /// overlap is returned as `Ok(None)`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fastq as fastq;
    /// use fastq::{Reader, Record};
    /// use libpairassembly::io::noodles::merge_pairs;
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
    /// // Merge the pairs into consensus reads.
    /// for result in merge_pairs(read_pairs) {
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

        let overlap_settings = OverlapParams::default();
        let validator = OverlapValidator::default();

        // Use the top-level checked assembler path.
        let assembler = Assembler::builder()
            .with_overlap_params(overlap_settings)
            .with_validator(validator)
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
    fn test_merge_pairs_with_perfect_overlap() {
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
    fn test_merge_pairs_mismatched_ids_yields_error() {
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
