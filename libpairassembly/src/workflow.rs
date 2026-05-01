pub use crate::prelude::*;

mod internal_demo {

    use crate::prelude::*;
    /// Demo for internal use showing how to run the `libpairassembly` API on a single pair of reads.
    /// The whole point of the library is to wrap it within your own iterator or stream-based system
    /// of progressively processing sequence reads, where merging occurs on each read pair.
    ///
    /// As such, the `libpairassembly` API is basically object-oriented in that routines are bundled
    /// with and performed per-object, and not on collections of objects or structs of arrays. That
    /// said, operations on each object could be seen as more data-oriented, as what is a sequence
    /// read object if not a struct of arrays of bytes?
    #[allow(dead_code, unused_variables)]
    fn demo() -> Result<()> {
        // Initialize some dummy reads
        let dummy_fwd_mate = SequenceRead::try_new("test", "ATGCC", "+4!;;")?;
        let dummy_rev_mate = SequenceRead::try_new("test", "AACTG", "+4!:;")?;

        // Pair up the reads. This will implicitly check that the two reads can be paired.
        let pair_input = PairInput::new(dummy_fwd_mate, dummy_rev_mate);

        // Initialize settings for overlapping and for validating those overlaps. We'll just use
        // defaults for demonstration purposes. Note that these are currently consumed, though this
        // may change in the future.
        let overlap_settings = OverlapParams::default();
        let validator = BaseCallValidator::default();

        // TODO: there should also be merge settings and correction settings, right?

        // The method chain below is it--the whole frontend of the library. `libpairassembly` presents an
        // API to be used as a chain of composable methods to execute the whole pair assembly pipeline.
        // It centers around four tasks that are mostly, though not always, separate concerns:
        // overlapping, validating overlaps, merging, and correction.
        //
        // Each step uses information from the last to improve and refine state before it's passed
        // onto the next step:
        //  - `.try_find_overlap()` finds the bounds of the best-available overlap between paired read, if any.
        //  - `.try_validate()` takes those bounds and interrogates it along with the reads that contain it to ensure that the overlap is good enough; if it is, it passes on the overlap bounds together with the underlying read data.
        //  - `.merge()` "interpolates" the paired reads around their validated overlap, but holds onto the original quality scores.
        //  - `.correct_quality_scores()` uses the original quality scores to adjust the merged quality scores to reflect whether the two reads agreed or disagreed about each base-call in the overlap. Agreeing base-calls should mean a quality score bump, and the opposite should occur for disagreeing base-calls.
        let assembler = Assembler::builder()
            .overlap(overlap_settings)
            .validate(validator)
            .build()?;

        let corrected_merged_read_yay = assembler
            .on_pair(&pair_input)?
            .overlap()?
            .validate()?
            .merge()?
            .correct()?
            .into_owned_read()?;

        let id = corrected_merged_read_yay.id();
        let new_sequence = corrected_merged_read_yay.sequence();
        let new_qualities = corrected_merged_read_yay.quality_scores();

        unimplemented!()
    }
}
