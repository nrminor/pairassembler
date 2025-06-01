use crate::{Result, SequenceRead, validate::ValidatedOverlap};

// pub trait Merge<'read> {
//     fn merge(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>>;
//     fn call_consensus_seq(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>>;
//     fn fwd_source_qual(&self) -> &[u8];
//     fn rev_source_qual(&self) -> &[u8];
// }

impl<'read> ValidatedOverlap<'read> {
    /// A core, heavy-lifter function in this crate. `merge()` takes views into the original reads
    /// as well as their overlaps and generates a consensus read. It notably does not, on its own,
    /// handle quality score corrections yet. Instead, its main role is to index into the forward
    /// and reverse read mates, concatenate references to their bytes into contiguous slices, and
    /// decide which base-call should be selected for each base. Given that, this method could in
    /// principle be used for sequence read formats without quality scores like FASTA. It also allows
    /// users to implement their own merging logic to improve performance or add other enhancements
    /// independent of how correction is implemented.
    pub fn merge(&self) -> Result<UncorrectedMergedRead<'read>> {
        // TODO: this function is much too long
        self.check_lengths()?;

        // TODO: The following section could be moved to its own default implementation trait once
        // the necessary methods for an upstream `Overlap` trait is implemented
        // -----------------------------------------------------------------------------------------
        let old_fwd = &self.mates.fwd_mate;
        let old_rev = &self.mates.rev_mate;

        // pull out a view to the forward overhang
        let LeftOverhang {
            start: left_overhang_start,
            stop: left_overhang_end,
            seq: left_overhang_seq,
            qual: left_overhang_qual,
        } = self.get_left_overhang(old_fwd, old_rev);

        // do a rotation on the reverse bounds to find the bases that belong on the reverse overhang
        let RightOverhang {
            start: right_overhang_start,
            stop: right_overhang_end,
            seq: right_overhang_seq,
            qual: right_overhang_qual,
        } = self.get_right_overhang(old_fwd, old_rev);

        // rebundle bases and quality scores to help choose which base should be at each position in
        // the consensus. Separate references to bytes in slices will become iterators of tuple pairs
        let fwd_entrants = self
            .overlap
            .r1_seq_view
            .iter()
            .zip(self.overlap.r1_qual_view.iter());
        let rev_entrants = self
            .overlap
            .r2_seq_view
            .iter()
            .zip(self.overlap.r2_qual_view.iter());
        // -----------------------------------------------------------------------------------------

        // choose the best base at each position in the overlap and fill it into an fixed size array
        let (consensus_seq, consensus_qual): (Vec<u8>, Vec<u8>) = fwd_entrants
            .zip(rev_entrants)
            .map(|((fwd_base, fwd_qual), (rev_base, rev_qual))| {
                match (fwd_base, fwd_qual, rev_base, rev_qual) {
                    _ if fwd_qual > rev_qual => (*fwd_base, *fwd_qual),
                    _ if fwd_qual < rev_qual => (*rev_base, *rev_qual),

                    // TODO: How best to break ties?
                    //
                    // This is a surprisingly tricky question to answer. Some methods do what we do here:
                    // just choose the forward mate's base, which is simple and deterministic, but can
                    // also introduce a consistent bias toward base-calls in the forward read. The same
                    // would be true if you arbitrarily used the reverse mate, and in either case, if
                    // you're processing Illumina reads, you're likely choosing the base that's in the
                    // lower-quality tail toward the end of each read.
                    //
                    // Another approach is to choose one or the other randomly, which does away with
                    // artificial bias, but also does away with determinism (I know I know, "just use
                    // a seed!" You could and should, but unfortunately, not to burst your bubble, but
                    // seeds do not always guarantee determinism, as they can behave differently on
                    // different platforms, OS's, etc. For example, 50-year old C code that set an
                    // explicit seed before pseudo-random number generation probably started at a
                    // different number when it was written than it does not). This feels dirty to me,
                    // and given the likely marginal frequency of these ties, I've decided against it.
                    // That said, if the world moves in the direction of reduced-bin quality scores
                    // (e.g., only having two quality scores corresponding to pass or fail so qualities
                    // can be encoded with two bits, or something similar), the random choice may
                    // become more attractive because ties will become the normal case and not the
                    // edge case.
                    //
                    // A third approach is to use available information from the basecalls prior to
                    // the base-call currently under examination, but this is a can of worms too.
                    // Illumina reads famously decline in quality from where the polymerase starts to
                    // when it falls of the strand, with a sharp drop in quality toward the end. Paired
                    // reading of the same short DNA strand usually means that overlaps between the two
                    // are at the ends of the forward mate, the reverse mate, or both. Depending on how
                    // sharp the increase in terminal low-quality is, it could be quite reasonable to
                    // distrust a call the more bases of high-quality precede it.
                    //
                    // You could of course innovate your way out of this dilemma, as there's still
                    // available information to use. For example, you could only look at the average
                    // quality of bases in a short window preceding the current base-call, e.g., 4 or 5
                    // calls back. This local average might throw away less information than a gross
                    // mean of the whole read so far (which is an example of the general phenomenon of
                    // summary statistics throwing away enough information on variance or data
                    // structure as to be useless for large, heterogeneous datasets, but I digress...).
                    //
                    // You could also adjust average quality by the number of bases assessed such that
                    // the trust in the average decreases non-linearly as the number of bases examined
                    // goes up. This may be the best of all worlds. Someone should implement that! 😉
                    _ if fwd_qual == rev_qual => (*fwd_base, *fwd_qual),

                    // hey man someone should do something here
                    _ => unimplemented!(),
                }
            })
            .collect();

        // make the overhang - consensus overlap - overhang sandwich, while also performing a sanity
        // check on overhangs and merged output if in debug testing
        let (full_consensus_seq, full_consensus_qual) = utils::concat_full_consensus(
            left_overhang_seq,
            &consensus_seq,
            &right_overhang_seq,
            left_overhang_qual,
            &consensus_qual,
            &right_overhang_qual,
        )?;

        let uncorrected_merged_read = UncorrectedMergedRead {
            // TODO: Handling ID's and full headers should be more sophisticated. In its current form
            // it will be lossy.
            id: old_fwd.id().to_owned(),
            consensus_seq: full_consensus_seq,
            consensus_qual: full_consensus_qual,
            fwd_source_seq: self.overlap.r1_seq_view,
            fwd_source_qual: self.overlap.r1_qual_view,
            // TODO: Unfortunately, these have to be owned and therefore must be cloned. Future
            // optimizations may be able to remove the clones. Ultimately, the the issue is that
            // we've moved ownership of each base and quality into the new consensus, meaning that
            // if we want to also look at the old values for correction, we'll need clones.
            rev_source_seq: self.overlap.r2_seq_view.clone(),
            rev_source_qual: self.overlap.r2_qual_view.clone(),
        };

        Ok(uncorrected_merged_read)
    }

    fn get_left_overhang<'overlap, 'mate>(
        &'overlap self,
        old_fwd: &'mate SequenceRead<'mate>,
        old_rev: &'mate SequenceRead<'mate>,
    ) -> LeftOverhang<'mate>
    where
        'overlap: 'mate,
    {
        // pull out a view to the forward overhang
        let start = 0;
        let stop = self.overlap.r1_start_offset;
        let seq = &old_fwd.sequence().as_bytes()[start..stop];
        let qual = &old_fwd.quality_scores().as_bytes()[start..stop];

        LeftOverhang {
            start,
            stop,
            seq,
            qual,
        }
    }

    fn get_right_overhang<'overlap, 'mate>(
        &'overlap self,
        old_fwd: &'mate SequenceRead<'mate>,
        old_rev: &'mate SequenceRead<'mate>,
    ) -> RightOverhang
    where
        'overlap: 'mate,
    {
        let start = self.overlap.r2_start_offset;
        let stop = old_rev.len() - 1;
        let seq = old_rev
            .reverse_complement()
            .drain(start..stop)
            .collect::<Vec<_>>();
        let qual = {
            // let's do some temporary in-place mutation to limit allocating the same information in
            // multiple Vec's all over the heap
            let mut initial_vec = old_rev.quality_scores().as_bytes().to_vec();
            initial_vec.reverse();
            initial_vec.drain(start..stop);
            initial_vec
        };

        RightOverhang {
            start,
            stop,
            seq,
            qual,
        }
    }

    fn call_consensus_seq(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>> {
        unimplemented!()
    }

    fn fwd_source_qual(&self) -> &[u8] {
        unimplemented!()
    }

    fn rev_source_qual(&self) -> &[u8] {
        unimplemented!()
    }
}

#[derive(Debug)]
pub struct UncorrectedMergedRead<'overlap> {
    pub id: String,
    pub consensus_seq: Vec<u8>,
    pub consensus_qual: Vec<u8>,
    pub fwd_source_seq: &'overlap [u8],
    pub fwd_source_qual: &'overlap [u8],
    pub rev_source_seq: Vec<u8>,
    pub rev_source_qual: Vec<u8>,
}

impl UncorrectedMergedRead<'_> {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn sequence(&self) -> &[u8] {
        self.consensus_seq.as_slice()
    }

    pub fn sequence_owned(self) -> Vec<u8> {
        self.consensus_seq
    }

    pub fn qualities(&self) -> &[u8] {
        self.consensus_qual.as_slice()
    }

    pub fn qualities_owned(self) -> Vec<u8> {
        self.consensus_qual
    }
}

struct LeftOverhang<'a> {
    start: usize,
    stop: usize,
    seq: &'a [u8],
    qual: &'a [u8],
}

struct RightOverhang {
    start: usize,
    stop: usize,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

// use utils::*;
mod utils {
    use crate::{
        Result,
        errors::MergeError::{MergedLengthMismatch, MismatchedConsensusSliceLengths},
        validate::ValidatedOverlap,
    };

    impl ValidatedOverlap<'_> {
        pub(super) fn check_lengths(&self) -> Result<()> {
            // TODO: All these lengths should probably be computed once and stored in the struct
            let r1_seq_len = self.overlap.r1_seq_view.len();
            let r2_seq_len = self.overlap.r2_seq_view.len();
            let r1_qual_len = self.overlap.r1_qual_view.len();
            let r2_qual_len = self.overlap.r2_qual_view.len();

            let r1_seq_qual_match = r1_seq_len == r1_qual_len;
            debug_assert!(
                r1_seq_qual_match,
                "Length mismatch in r1 sequence and quality views"
            );
            let r2_seq_qual_match = r2_seq_len == r2_qual_len;
            debug_assert!(
                r2_seq_qual_match,
                "Length mismatch in r2 sequence and quality views"
            );

            let r1_r2_len_check = r1_seq_len == r2_seq_len;
            debug_assert!(r1_r2_len_check, "Overlap sequences are not the same length");
            if !r1_r2_len_check {
                return Err(MismatchedConsensusSliceLengths {
                    fwd_len: r1_seq_len,
                    rev_len: r2_seq_len,
                }
                .into());
            }

            Ok(())
        }
    }

    pub(super) fn concat_full_consensus(
        left_overhang_seq: &[u8],
        consensus_seq: &[u8],
        right_overhang_seq: &[u8],
        left_overhang_qual: &[u8],
        consensus_qual: &[u8],
        right_overhang_qual: &[u8],
    ) -> Result<(Vec<u8>, Vec<u8>)> {
        let post_merge_len =
            left_overhang_seq.len() + consensus_seq.len() + right_overhang_seq.len();

        let full_consensus_seq = {
            let mut seq_sandwich = Vec::with_capacity(post_merge_len);
            seq_sandwich.extend_from_slice(left_overhang_seq);
            seq_sandwich.extend_from_slice(consensus_seq);
            seq_sandwich.extend_from_slice(right_overhang_seq);
            seq_sandwich
        };
        let full_consensus_qual = {
            let mut qual_sandwich = Vec::with_capacity(post_merge_len);
            qual_sandwich.extend_from_slice(left_overhang_qual);
            qual_sandwich.extend_from_slice(consensus_qual);
            qual_sandwich.extend_from_slice(right_overhang_qual);
            qual_sandwich
        };
        debug_assert_eq!(full_consensus_seq.len(), full_consensus_qual.len());
        debug_assert!(full_consensus_seq.len() >= consensus_seq.len());

        let length_check = full_consensus_seq.len() == post_merge_len;
        debug_assert!(
            length_check,
            "Post-merge sequence length computation mismatch"
        );
        if !length_check {
            return Err(MergedLengthMismatch {
                expected: post_merge_len,
                actual: full_consensus_seq.len(),
            }
            .into());
        }

        Ok((full_consensus_seq, full_consensus_qual))
    }
}
