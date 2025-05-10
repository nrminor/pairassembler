use crate::{ReadMates, SequenceRead, ValidatedOverlap, overlap};

pub trait Merge<'read> {
    fn merge(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>>;
    fn call_consensus_seq(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>>;
    fn fwd_source_qual(&self) -> &[u8];
    fn rev_source_qual(&self) -> &[u8];
}

impl<'read> Merge<'read> for ValidatedOverlap<'read> {
    /// A core, heavy-lifter function in this crate. `merge()` takes views into the original reads
    /// as well as their overlaps and generates a consensus read. It notably does not, on its own,
    /// handle quality score corrections yet. Instead, its main role is to index into the forward
    /// and reverse read mates, concatenate references to their bytes into contiguous slices, and
    /// decide which base-call should be selected for each base. Given that, this method could in
    /// principle be used for sequence read formats without quality scores like FASTA. It also allows
    /// users to implement their own merging logic to improve performance or add other enhancements
    /// independent of how correction is implemented.
    fn merge(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>> {
        // TODO: The following section should be moved to its own default implementation trait once
        // the necessary methods for an upstream `Overlap` trait is implemented
        // -----------------------------------------------------------------------------------------
        let old_fwd = &self.mates.fwd_mate;
        let old_rev = &self.mates.rev_mate;

        // pull out a view to the forward overhang
        let fwd_overhang_start = 0;
        let fwd_overhang_end = self.overlap.r1_start_offset;
        let left_overhang = &old_fwd.sequence().as_bytes()[0..fwd_overhang_start];

        // do a rotation on the reverse bounds to find the bases that belong on the reverse overhang
        let right_overhang_start = self.overlap.r2_start_offset;
        let right_overhang_end = old_rev.len() - 1;
        let right_overhang =
            &old_rev.reverse_complement()[right_overhang_start..right_overhang_end];

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
        let consensus_seq = fwd_entrants
            .zip(rev_entrants)
            .map(|((fwd_base, fwd_qual), (rev_base, rev_qual))| {
                match (fwd_base, fwd_qual, rev_base, rev_qual) {
                    _ if fwd_qual > rev_qual => *fwd_base,
                    _ if fwd_qual < rev_qual => *rev_base,

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
                    _ if fwd_qual == rev_qual => *fwd_base,
                    _ => unimplemented!(),
                }
            })
            .collect::<Vec<_>>();
        let overlap = &consensus_seq;

        // make the overhang - consensus overlap - overhang sandwich
        let post_merge_len = left_overhang.len() + consensus_seq.len() + right_overhang.len();
        let seq_sandwich = {
            let mut seq_sandwich = Vec::with_capacity(post_merge_len);
            seq_sandwich.extend_from_slice(left_overhang);
            seq_sandwich.extend_from_slice(&consensus_seq);
            seq_sandwich.extend_from_slice(right_overhang);
            seq_sandwich
        };

        todo!()
    }

    fn call_consensus_seq(&self) -> color_eyre::Result<UncorrectedMergedRead<'read>> {
        todo!()
    }

    fn fwd_source_qual(&self) -> &[u8] {
        todo!()
    }

    fn rev_source_qual(&self) -> &[u8] {
        todo!()
    }
}

#[derive(Debug)]
pub(crate) struct UncorrectedMergedRead<'overlap> {
    pub(crate) id: String,
    pub(crate) consensus_seq: Vec<u8>,
    pub(crate) consensus_qual: Vec<u8>,
    pub(crate) fwd_source_seq: &'overlap [u8],
    pub(crate) fwd_source_qual: &'overlap [u8],
    pub(crate) rev_source_seq: &'overlap [u8],
    pub(crate) rev_source_qual: &'overlap [u8],
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

use helpers::*;
mod helpers {}
