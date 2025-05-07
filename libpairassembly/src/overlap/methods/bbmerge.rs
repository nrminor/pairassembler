#![allow(clippy::pedantic, clippy::perf)]
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

/*!
This module is a very stripped down refactor of Brian Bushnell's BBMerge paired-read
merging method. The point of this module is to extract the core merging logic that operates
per-read, including error correction. In so doing, it will cut away all the sophisticated
logging BBMerge does, which involves a large number of tickers and global (well, BBMerge-class-
level) boolean checks.

This module will also differ in that joining and error-correction of merged reads will not
be optional, nor will they be mutually exclusive as they are with BBMerge. It will also not
handle adapter sequences or trim so that only the overlap remains. All it will do is find
overlaps, if any, join overlapping bases, and perform error correction.

As such, the majority of the logic here falls within the if-statement in `BBMerge.java` that
tests whether a merge was detected: `if (bestInsert > 0) { ... }`, as well as within `.parseInsert()`
and `.processReadPath()` methods.

Eventually, this system, like the FASTP implementation, will be plugged into an asynchronous
streaming engine FASTQ records and IO, though this system may be implemented in the command-line
tool associated with this library, and not in the library itself.
*/

use color_eyre::Result;

#[derive(Debug)]
struct Read<'a> {
    id: &'a str,
    seq: &'a str,
    qual: &'a str,
}

impl<'a> Read<'a> {
    fn new(id: &'a str, seq: &'a str, qual: &'a str) -> Self {
        Self { id, seq, qual }
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.seq.len(), self.qual.len());
        self.seq.len()
    }

    pub fn reverse_complement(&self) -> String {
        todo!()
    }
}

#[derive(Debug)]
struct ReadMates<'read> {
    fwd_mate: Read<'read>,
    rev_mate: Read<'read>,
    insert_min: usize,
    insert_max: usize,
}

impl ReadMates<'_> {
    pub fn parse_insert(&self) -> Result<usize> {
        todo!()
    }

    pub async fn process_read_pair(&self) -> Result<usize> {
        todo!()
    }

    pub async fn join_reads(&self, insert: usize) -> Result<Read> {
        assert!(insert > 0);

        let length_sum = self.fwd_mate.len() + self.rev_mate.len();
        let overlap = insert.min(length_sum - insert);

        let mut mismatches = 0;
        let mut start = 0;
        let mut stop = 0;

        todo!()
    }

    pub async fn error_correct_with_insert(&self, insert: usize) -> Result<usize> {
        if insert >= self.fwd_mate.len() + self.rev_mate.len() {
            return Ok(0);
        }

        assert!(insert > 0);

        let joined = self.join_reads(insert).await?;

        assert_ne!(joined.len(), 0);

        let len_joined = joined.len();
        let limit1 = len_joined.min(self.fwd_mate.len());
        let limit2 = len_joined.min(self.rev_mate.len());

        let old_fwd_bases = self.fwd_mate.seq.as_bytes();
        let old_rev_bases = self.rev_mate.seq.as_bytes();

        let mut fwd_errors = 0;
        let mut rev_errors = 0;

        let new_fwd_bases = &joined.seq.as_bytes()[0..limit1];
        for i in 0..new_fwd_bases.len() {
            if old_fwd_bases[i] != new_fwd_bases[i] {
                fwd_errors += 1;
            }
        }

        let new_rev_bases = &joined.seq.as_bytes()[0..limit1];
        let loop_offset = old_rev_bases.len() - new_rev_bases.len();
        for (i, &base) in new_rev_bases.iter().enumerate() {
            let j = loop_offset + i;
            if old_rev_bases[j] != base {
                rev_errors += 1;
            }
        }

        // As I refactor this, it's looking like `ecco` doesn't actually tell
        // BBMerge to do any error-correction. It just tells it to tally up the
        // error's it has already corrected so they can be logged out.
        //
        // That said, the error-correction may be performed downstream closer
        // to writing out.

        Ok(fwd_errors + rev_errors)
    }

    pub async fn find_overlap(&self, ecco: bool) -> Result<Option<i32>> {
        let true_size = self.parse_insert()?;
        let initial_len1 = self.fwd_mate.len();
        let initial_len2 = self.rev_mate.len();
        let min_read_len = initial_len1.min(initial_len2);
        let len_sum = initial_len1 + initial_len2;
        let true_overlap = true_size.min(min_read_len).min(len_sum - true_size);

        let best_insert = self.process_read_pair().await?;
        // debug!(
        //     "The calculated true size for {:?} was {:?}, and the best insert was {:?}",
        //     true_size, self, best_insert
        // );

        if best_insert > 0 {
            // JAVA WE WON'T CONVERT BUT I'M HOLDING ONTO FOR NOW
            // if (trueOverlap > 0 && trueOverlap <= 8 && EXEMPT_SHORT_OVERLAPS_FROM_STATS) {
            // 	// Temporarily ignore for purpose of training
            // } else if (bestInsert == trueSize) {
            // 	correctCount++;
            // 	insertSumCorrect += bestInsert;
            // } else {
            // 	incorrectCount++;
            // 	insertSumIncorrect += bestInsert;
            // }

            // r1.setInsert(bestInsert);

            let insert_min = best_insert.min(self.insert_min);
            let insert_max = best_insert.min(self.insert_max);

            // JAVA WE WON'T CONVERT BUT I'M HOLDING ONTO FOR NOW
            // hist[Tools.min(bestInsert, hist.length - 1)]++;
            //
            // Bushnell also makes it an option NOT to join overlapping reads; here, we will
            // always join them

            let read2_revcomp = self.rev_mate.reverse_complement();

            let joined_read = self.join_reads(best_insert).await?;

            // Bushnell also gives error correcting as an option with the following java:
            //
            // ```
            // } else if (ecco) {
            // 	r2.reverseComplement();
            // 	errorCorrectWithInsert(r1, r2, bestInsert);
            // 	r2.reverseComplement();
            // 	errorsCorrectedT += (r1.errors + r2.errors);
            //
            // ```
            //
            // As before, error-correction will not be an option in this library; it
            // will always be performed. As such, the joing logic here will be more
            // akin to the BBMerge method called `.errorCorrectWithInsert()` than
            // `.joinReads()` alone.

            // There are also lots of counters getting incremented for runtime
            // introspection in BBMerge, like what follows. We will not use any
            // of these counters for not, though they're of course not a bad idea.
            //
            // ```
            // if (bestInsert == RET_AMBIG) {
            // 	ambiguousCount++;
            // } else if (bestInsert == RET_SHORT) {
            // 	tooShortCount++;
            // } else if (bestInsert == RET_LONG) {
            // 	tooLongCount++;
            // } else if (bestInsert == RET_NO_SOLUTION) {
            // 	noSolutionCount++;
            // }

            // if (trueSize < initialLen1 + initialLen2 && trueSize > 0) {
            // 	unmergedOverlappingCount++;
            // }

            // ```
        }

        todo!()
    }
}

#[derive(Debug)]
pub struct BBMerge {
    method: MergeMethod,
    ecco: bool,
}

#[derive(Debug, Default)]
pub enum MergeMethod {
    Mapping,
    #[default]
    MateByOverlap,
}

impl BBMerge {
    pub async fn run(&self) -> Result<()> {
        self.process_reads().await
    }

    pub async fn process_reads(&self) -> Result<()> {
        todo!()
    }
}
