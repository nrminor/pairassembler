use color_eyre::eyre::eyre;

pub use crate::prelude::*;

#[allow(dead_code, unused_variables)]
pub fn test() -> color_eyre::Result<()> {
    // Initialize some dummy reads
    let dummy_fwd_mate = SequenceRead::try_new("test", "ATGCC", "+4!;;")?;
    let dummy_rev_mate = SequenceRead::try_new("test", "AACTG", "+4!:;")?;

    // Pair up the reads. This will implicitly check that the two reads can be paired.
    let mates = ReadMates::from(dummy_fwd_mate, dummy_rev_mate)?;

    // Initialize settings for overlapping and for validating those overlaps. We'll just use
    // defaults for demonstration purposes. Note that these are currently consumed, though this
    // may change in the future.
    let overlap_settings = OverlapParams::default();
    let validator = BaseCallValidator::default();

    // Use a chain of methods to execute the whole pipeline
    let corrected_merged_read_yay = mates
        .try_find_overlap(&overlap_settings)?
        .ok_or_else(|| eyre!("Overlapping failed."))?
        .try_validate(&mates, &validator)?
        .merge()?
        .correct_quality_scores()?;

    let id = corrected_merged_read_yay.id();
    let new_sequence = corrected_merged_read_yay.sequence();
    let new_qualities = corrected_merged_read_yay.qualities();

    unimplemented!()
}
