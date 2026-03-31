use crate::{
    Result, SequenceRead,
    assembler::HasMergeableOverlap,
    correct::CorrectedMergedRead,
    errors::MergeError::{MergedLengthMismatch, MismatchedConsensusSliceLengths},
    validate::ValidatedOverlap,
};

#[derive(Debug, Clone)]
pub struct MergedRead {
    id: String,
    consensus_seq: Vec<u8>,
    consensus_qual: Vec<u8>,
    provenance: MergeProvenance,
}

#[derive(Debug, Clone)]
pub struct MergeProvenance {
    overlap_len: usize,
    fwd_overlap_seq: Vec<u8>,
    fwd_overlap_qual: Vec<u8>,
    rev_overlap_seq: Vec<u8>,
    rev_overlap_qual: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct MergeView {
    pub id: String,
    pub left_seq: Vec<u8>,
    pub left_qual: Vec<u8>,
    pub fwd_overlap_seq: Vec<u8>,
    pub fwd_overlap_qual: Vec<u8>,
    pub rev_overlap_seq: Vec<u8>,
    pub rev_overlap_qual: Vec<u8>,
    pub right_seq: Vec<u8>,
    pub right_qual: Vec<u8>,
}

impl MergedRead {
    pub(crate) fn from_parts(
        id: String,
        consensus_seq: Vec<u8>,
        consensus_qual: Vec<u8>,
        provenance: MergeProvenance,
    ) -> Self {
        Self {
            id,
            consensus_seq,
            consensus_qual,
            provenance,
        }
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn sequence(&self) -> &[u8] {
        self.consensus_seq.as_slice()
    }

    pub fn qualities(&self) -> &[u8] {
        self.consensus_qual.as_slice()
    }

    pub fn len(&self) -> usize {
        self.consensus_seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.consensus_seq.is_empty()
    }

    pub fn sequence_owned(self) -> Vec<u8> {
        self.consensus_seq
    }

    pub fn qualities_owned(self) -> Vec<u8> {
        self.consensus_qual
    }

    pub fn provenance(&self) -> &MergeProvenance {
        &self.provenance
    }

    pub fn correct(self) -> Result<CorrectedMergedRead> {
        self.into_uncorrected().correct()
    }

    pub(crate) fn into_uncorrected(self) -> UncorrectedMergedRead {
        let MergedRead {
            id,
            consensus_seq,
            consensus_qual,
            provenance,
        } = self;

        UncorrectedMergedRead::from_parts(
            id,
            consensus_seq,
            consensus_qual,
            provenance.fwd_overlap_seq,
            provenance.fwd_overlap_qual,
            provenance.rev_overlap_seq,
            provenance.rev_overlap_qual,
        )
    }
}

impl MergeProvenance {
    pub(crate) fn from_parts(
        overlap_len: usize,
        fwd_overlap_seq: Vec<u8>,
        fwd_overlap_qual: Vec<u8>,
        rev_overlap_seq: Vec<u8>,
        rev_overlap_qual: Vec<u8>,
    ) -> Self {
        Self {
            overlap_len,
            fwd_overlap_seq,
            fwd_overlap_qual,
            rev_overlap_seq,
            rev_overlap_qual,
        }
    }

    pub fn overlap_len(&self) -> usize {
        self.overlap_len
    }

    pub fn fwd_overlap_seq(&self) -> &[u8] {
        self.fwd_overlap_seq.as_slice()
    }

    pub fn fwd_overlap_qual(&self) -> &[u8] {
        self.fwd_overlap_qual.as_slice()
    }

    pub fn rev_overlap_seq(&self) -> &[u8] {
        self.rev_overlap_seq.as_slice()
    }

    pub fn rev_overlap_qual(&self) -> &[u8] {
        self.rev_overlap_qual.as_slice()
    }
}

impl crate::assembler::IntoOwnedRecordParts for MergedRead {
    fn into_owned_record_parts(self) -> (String, Vec<u8>, Vec<u8>) {
        (self.id, self.consensus_seq, self.consensus_qual)
    }
}

pub fn merge_from<T>(input: &T) -> Result<MergedRead>
where
    T: HasMergeableOverlap,
{
    let view = input.merge_view()?;
    merge_kernel(view)
}

fn merge_kernel(view: MergeView) -> Result<MergedRead> {
    let overlap_len = view.fwd_overlap_seq.len();
    if overlap_len != view.rev_overlap_seq.len() {
        return Err(MismatchedConsensusSliceLengths {
            fwd_len: overlap_len,
            rev_len: view.rev_overlap_seq.len(),
        }
        .into());
    }

    let mut consensus_overlap_seq = Vec::with_capacity(overlap_len);
    let mut consensus_overlap_qual = Vec::with_capacity(overlap_len);

    for i in 0..overlap_len {
        let fb = view.fwd_overlap_seq[i];
        let fq = view.fwd_overlap_qual[i];
        let rb = view.rev_overlap_seq[i];
        let rq = view.rev_overlap_qual[i];

        if fq >= rq {
            consensus_overlap_seq.push(fb);
            consensus_overlap_qual.push(fq);
        } else {
            consensus_overlap_seq.push(rb);
            consensus_overlap_qual.push(rq);
        }
    }

    let expected_len = view.left_seq.len() + overlap_len + view.right_seq.len();
    let mut full_seq = Vec::with_capacity(expected_len);
    full_seq.extend_from_slice(&view.left_seq);
    full_seq.extend_from_slice(&consensus_overlap_seq);
    full_seq.extend_from_slice(&view.right_seq);

    let mut full_qual = Vec::with_capacity(expected_len);
    full_qual.extend_from_slice(&view.left_qual);
    full_qual.extend_from_slice(&consensus_overlap_qual);
    full_qual.extend_from_slice(&view.right_qual);

    if full_seq.len() != expected_len {
        return Err(MergedLengthMismatch {
            expected: expected_len,
            actual: full_seq.len(),
        }
        .into());
    }

    let provenance = MergeProvenance::from_parts(
        overlap_len,
        view.fwd_overlap_seq,
        view.fwd_overlap_qual,
        view.rev_overlap_seq,
        view.rev_overlap_qual,
    );

    Ok(MergedRead::from_parts(
        view.id, full_seq, full_qual, provenance,
    ))
}

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
    ///
    /// # Errors
    ///
    /// Returns an error when overlap-derived bounds are inconsistent with mate sequence/quality
    /// lengths.
    pub fn merge(&self) -> Result<UncorrectedMergedRead> {
        // TODO: this function is much too long
        self.check_lengths()?;

        // TODO: The following section could be moved to its own default implementation trait once
        // the necessary methods for an upstream `Overlap` trait is implemented
        // -----------------------------------------------------------------------------------------
        let old_fwd = &self.read_pair().fwd_mate;
        let old_rev = &self.read_pair().rev_mate;

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
        let overlap = self.overlap();
        let fwd_entrants = overlap
            .forward_sequence()
            .iter()
            .zip(overlap.forward_qualities().iter());
        let rev_entrants = overlap
            .reverse_sequence()
            .iter()
            .zip(overlap.reverse_qualities().iter());
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

        let uncorrected_merged_read = UncorrectedMergedRead::from_parts(
            // TODO: Handling ID's and full headers should be more sophisticated. In its current form
            // it will be lossy.
            old_fwd.id().to_owned(),
            full_consensus_seq,
            full_consensus_qual,
            overlap.forward_sequence().to_vec(),
            overlap.forward_qualities().to_vec(),
            // TODO: Unfortunately, these have to be owned and therefore must be cloned. Future
            // optimizations may be able to remove the clones. Ultimately, the the issue is that
            // we've moved ownership of each base and quality into the new consensus, meaning that
            // if we want to also look at the old values for correction, we'll need clones.
            overlap.reverse_sequence().to_vec(),
            overlap.reverse_qualities().to_vec(),
        );

        Ok(uncorrected_merged_read)
    }

    fn get_left_overhang<'overlap, 'mate>(
        &'overlap self,
        old_fwd: &'mate SequenceRead<'mate>,
        old_rev: &'mate SequenceRead<'mate>,
    ) -> LeftOverhang<'mate>
    where
        'mate: 'overlap,
    {
        // pull out a view to the forward overhang
        let start = 0;
        let stop = self.overlap().forward_start_offset();
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
        'mate: 'overlap,
    {
        let start = self.overlap().reverse_end_offset() + 1;
        let stop = old_rev.len();
        let seq = old_rev
            .reverse_complement()
            .drain(start..stop)
            .collect::<Vec<_>>();
        let qual = {
            // let's do some temporary in-place mutation to limit allocating the same information in
            // multiple Vec's all over the heap
            let mut initial_vec = old_rev.quality_scores().as_bytes().to_vec();
            initial_vec.reverse();
            initial_vec.drain(start..stop).collect::<Vec<_>>()
        };

        RightOverhang {
            start,
            stop,
            seq,
            qual,
        }
    }

    fn call_consensus_seq(&self) -> color_eyre::Result<UncorrectedMergedRead> {
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
pub struct UncorrectedMergedRead {
    id: String,
    consensus_seq: Vec<u8>,
    consensus_qual: Vec<u8>,
    fwd_source_seq: Vec<u8>,
    fwd_source_qual: Vec<u8>,
    rev_source_seq: Vec<u8>,
    rev_source_qual: Vec<u8>,
}

impl UncorrectedMergedRead {
    pub(crate) fn from_parts(
        id: String,
        consensus_seq: Vec<u8>,
        consensus_qual: Vec<u8>,
        fwd_source_seq: Vec<u8>,
        fwd_source_qual: Vec<u8>,
        rev_source_seq: Vec<u8>,
        rev_source_qual: Vec<u8>,
    ) -> Self {
        Self {
            id,
            consensus_seq,
            consensus_qual,
            fwd_source_seq,
            fwd_source_qual,
            rev_source_seq,
            rev_source_qual,
        }
    }

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

    pub fn forward_source_seq(&self) -> &[u8] {
        self.fwd_source_seq.as_slice()
    }

    pub fn forward_source_qual(&self) -> &[u8] {
        self.fwd_source_qual.as_slice()
    }

    pub fn reverse_source_seq(&self) -> &[u8] {
        self.rev_source_seq.as_slice()
    }

    pub fn reverse_source_qual(&self) -> &[u8] {
        self.rev_source_qual.as_slice()
    }

    pub(crate) fn into_correction_parts(
        self,
    ) -> (String, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
        (
            self.id,
            self.consensus_seq,
            self.fwd_source_seq,
            self.fwd_source_qual,
            self.rev_source_seq,
            self.rev_source_qual,
        )
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
            let overlap = self.overlap();
            let r1_seq_len = overlap.forward_sequence().len();
            let r2_seq_len = overlap.reverse_sequence().len();
            let r1_qual_len = overlap.forward_qualities().len();
            let r2_qual_len = overlap.reverse_qualities().len();

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

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::merge_from;
    use crate::{PairOverlap, ReadPair, SequenceRead, validate::ValidatedOverlap};
    use proptest::{collection::vec, prelude::*};

    fn dna_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        vec(
            prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
            min_len..=max_len,
        )
        .prop_map(|chars| chars.into_iter().collect())
    }

    fn qual_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
        vec(33u8..=73u8, min_len..=max_len)
            .prop_map(|bytes| bytes.into_iter().map(char::from).collect())
    }

    #[derive(Debug, Clone)]
    struct MergeFixture {
        left_seq: String,
        overlap_fwd_seq: String,
        overlap_rev_seq: String,
        right_seq: String,
        left_qual: String,
        overlap_fwd_qual: String,
        overlap_rev_qual: String,
        right_qual: String,
    }

    prop_compose! {
        fn merge_fixture_strategy()
            (left_len in 0usize..=16, overlap_len in 4usize..=24, right_len in 0usize..=16)
            (
                left_seq in dna_string_strategy(left_len, left_len),
                overlap_fwd_seq in dna_string_strategy(overlap_len, overlap_len),
                overlap_rev_seq in dna_string_strategy(overlap_len, overlap_len),
                right_seq in dna_string_strategy(right_len, right_len),
                left_qual in qual_string_strategy(left_len, left_len),
                overlap_fwd_qual in qual_string_strategy(overlap_len, overlap_len),
                overlap_rev_qual in qual_string_strategy(overlap_len, overlap_len),
                right_qual in qual_string_strategy(right_len, right_len),
            ) -> MergeFixture
        {
            MergeFixture {
                left_seq,
                overlap_fwd_seq,
                overlap_rev_seq,
                right_seq,
                left_qual,
                overlap_fwd_qual,
                overlap_rev_qual,
                right_qual,
            }
        }
    }

    fn reverse_complement_dna(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|base| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                invalid => panic!("invalid DNA base in merge test fixture: {invalid}"),
            })
            .collect()
    }

    fn build_validated_overlap_from_rev_rc_parts(
        left_seq: &str,
        overlap_fwd_seq: &str,
        overlap_rev_seq: &str,
        right_seq: &str,
        left_qual: &str,
        overlap_fwd_qual: &str,
        overlap_rev_qual: &str,
        right_qual: &str,
    ) -> ValidatedOverlap<'static> {
        let fwd_seq = format!("{left_seq}{overlap_fwd_seq}");
        let fwd_qual = format!("{left_qual}{overlap_fwd_qual}");

        let rev_rc_seq = format!("{overlap_rev_seq}{right_seq}");
        let rev_rc_qual = format!("{overlap_rev_qual}{right_qual}");

        let rev_seq = reverse_complement_dna(&rev_rc_seq);
        let rev_qual = rev_rc_qual.chars().rev().collect::<String>();

        let fwd_static: &'static str = Box::leak(fwd_seq.into_boxed_str());
        let fwd_qual_static: &'static str = Box::leak(fwd_qual.into_boxed_str());
        let rev_static: &'static str = Box::leak(rev_seq.into_boxed_str());
        let rev_qual_static: &'static str = Box::leak(rev_qual.into_boxed_str());

        let mates = ReadPair::from(
            SequenceRead::new("read1", fwd_static, fwd_qual_static),
            SequenceRead::new("read1", rev_static, rev_qual_static),
        )
        .expect("test fixtures should produce valid paired reads");
        let mates_ref: &'static ReadPair<'static> = Box::leak(Box::new(mates));

        let overlap_len = overlap_fwd_seq.len();
        let fwd_start = left_seq.len();
        let fwd_end = fwd_start + overlap_len - 1;
        let overlap = PairOverlap::from_components(
            overlap_len,
            fwd_start,
            fwd_end,
            0,
            overlap_len - 1,
            &mates_ref.fwd_sequence_bytes()[fwd_start..=fwd_end],
            &mates_ref.fwd_quality_bytes()[fwd_start..=fwd_end],
            overlap_rev_seq.as_bytes().to_vec(),
            overlap_rev_qual.as_bytes().to_vec(),
        );

        ValidatedOverlap::from_parts(mates_ref, overlap)
    }

    fn oracle_merge(
        left_seq: &str,
        overlap_fwd_seq: &str,
        overlap_rev_seq: &str,
        right_seq: &str,
        left_qual: &str,
        overlap_fwd_qual: &str,
        overlap_rev_qual: &str,
        right_qual: &str,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut overlap_seq = Vec::with_capacity(overlap_fwd_seq.len());
        let mut overlap_qual = Vec::with_capacity(overlap_fwd_qual.len());

        for (((fwd_base, fwd_q), rev_base), rev_q) in overlap_fwd_seq
            .bytes()
            .zip(overlap_fwd_qual.bytes())
            .zip(overlap_rev_seq.bytes())
            .zip(overlap_rev_qual.bytes())
        {
            if fwd_q >= rev_q {
                overlap_seq.push(fwd_base);
                overlap_qual.push(fwd_q);
            } else {
                overlap_seq.push(rev_base);
                overlap_qual.push(rev_q);
            }
        }

        let mut full_seq = Vec::new();
        full_seq.extend_from_slice(left_seq.as_bytes());
        full_seq.extend_from_slice(&overlap_seq);
        full_seq.extend_from_slice(right_seq.as_bytes());

        let mut full_qual = Vec::new();
        full_qual.extend_from_slice(left_qual.as_bytes());
        full_qual.extend_from_slice(&overlap_qual);
        full_qual.extend_from_slice(right_qual.as_bytes());

        (full_seq, full_qual)
    }

    #[test]
    fn test_merge_perfect_full_overlap_roundtrip() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = PairOverlap::from_components(
            5,
            4,
            8,
            0,
            4,
            &mates.fwd_sequence_bytes()[4..=8],
            &mates.fwd_quality_bytes()[4..=8],
            mates.rev_mate.reverse_complement(),
            mates.rev_mate.quality_scores().as_bytes().to_vec(),
        );
        let validated = ValidatedOverlap::from_parts(&mates, overlap);

        let merged = validated
            .merge()
            .expect("validated overlap should merge without bounds errors");

        assert_eq!(merged.id(), "read1");
        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.qualities(), b"IIIIIIIII");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_with_left_overhang_preserves_prefix() {
        let r1 = SequenceRead::new("read1", "TTTTACGTA", "IIIIIIIII");
        let r2 = SequenceRead::new("read1", "TACGT", "IIIII");
        let mates = ReadPair::from(r1, r2).expect("test fixture reads should share the same id");

        let overlap = PairOverlap::from_components(
            5,
            4,
            8,
            0,
            4,
            &mates.fwd_sequence_bytes()[4..=8],
            &mates.fwd_quality_bytes()[4..=8],
            mates.rev_mate.reverse_complement(),
            mates.rev_mate.quality_scores().as_bytes().to_vec(),
        );
        let validated = ValidatedOverlap::from_parts(&mates, overlap);

        let merged = validated
            .merge()
            .expect("validated overlap should merge without bounds errors");

        assert_eq!(merged.sequence(), b"TTTTACGTA");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_with_right_overhang_preserves_suffix() {
        let validated = build_validated_overlap_from_rev_rc_parts(
            "TT", "ACGT", "ACGT", "GG", "II", "IIII", "IIII", "II",
        );

        let merged = validated
            .merge()
            .expect("validated overlap with right overhang should merge");

        assert_eq!(merged.sequence(), b"TTACGTGG");
        assert_eq!(merged.qualities(), b"IIIIIIII");
        assert_eq!(merged.sequence().len(), merged.qualities().len());
    }

    #[test]
    fn test_merge_tie_on_quality_prefers_forward_base() {
        let validated = build_validated_overlap_from_rev_rc_parts(
            "", "AAAA", "TTTT", "", "", "IIII", "IIII", "",
        );

        let merged = validated
            .merge()
            .expect("equal-quality overlap should merge deterministically");

        assert_eq!(merged.sequence(), b"AAAA");
        assert_eq!(merged.qualities(), b"IIII");
    }

    #[test]
    fn test_merge_from_matches_legacy_validated_overlap_merge() {
        let validated = build_validated_overlap_from_rev_rc_parts(
            "TT", "ACGT", "TCGT", "GG", "II", "IIII", "IIII", "II",
        );

        let legacy = validated
            .merge()
            .expect("legacy validated-overlap merge should succeed");
        let migrated = merge_from(&validated)
            .expect("generic merge_from should match validated-overlap merge behavior");

        assert_eq!(migrated.id(), legacy.id());
        assert_eq!(migrated.sequence(), legacy.sequence());
        assert_eq!(migrated.qualities(), legacy.qualities());
        assert_eq!(
            migrated.provenance().fwd_overlap_seq(),
            legacy.forward_source_seq()
        );
        assert_eq!(
            migrated.provenance().fwd_overlap_qual(),
            legacy.forward_source_qual()
        );
        assert_eq!(
            migrated.provenance().rev_overlap_seq(),
            legacy.reverse_source_seq()
        );
        assert_eq!(
            migrated.provenance().rev_overlap_qual(),
            legacy.reverse_source_qual()
        );
    }

    proptest! {
        #[test]
        fn proptest_merge_matches_oracle_for_constructed_overlap(
            fixture in merge_fixture_strategy(),
        ) {
            let validated = build_validated_overlap_from_rev_rc_parts(
                &fixture.left_seq,
                &fixture.overlap_fwd_seq,
                &fixture.overlap_rev_seq,
                &fixture.right_seq,
                &fixture.left_qual,
                &fixture.overlap_fwd_qual,
                &fixture.overlap_rev_qual,
                &fixture.right_qual,
            );

            let merged = validated
                .merge()
                .expect("constructed overlap fixture should merge");
            let (expected_seq, expected_qual) = oracle_merge(
                &fixture.left_seq,
                &fixture.overlap_fwd_seq,
                &fixture.overlap_rev_seq,
                &fixture.right_seq,
                &fixture.left_qual,
                &fixture.overlap_fwd_qual,
                &fixture.overlap_rev_qual,
                &fixture.right_qual,
            );

            prop_assert_eq!(merged.sequence(), expected_seq.as_slice());
            prop_assert_eq!(merged.qualities(), expected_qual.as_slice());
            prop_assert_eq!(merged.sequence().len(), merged.qualities().len());
            prop_assert_eq!(merged.forward_source_seq(), fixture.overlap_fwd_seq.as_bytes());
            prop_assert_eq!(merged.forward_source_qual(), fixture.overlap_fwd_qual.as_bytes());
            prop_assert_eq!(merged.reverse_source_seq(), fixture.overlap_rev_seq.as_bytes());
            prop_assert_eq!(merged.reverse_source_qual(), fixture.overlap_rev_qual.as_bytes());
        }
    }
}
