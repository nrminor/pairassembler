use proptest::prelude::*;

use crate::{Assembler, OverlapParams, PairInput, SequenceRead, assembler::OverlapOutcome};

proptest! {
    #[test]
    fn proptest_paired_end_insert_geometry_discovers_expected_overlap(
        left_flank_len in 1usize..=48,
        overlap_len in 20usize..=96,
        right_flank_len in 1usize..=48,
    ) {
        let left_flank = "T".repeat(left_flank_len);
        let overlap = "A".repeat(overlap_len);
        let right_flank = "C".repeat(right_flank_len);

        let r1_sequence = format!("{left_flank}{overlap}");
        let oriented_r2_sequence = format!("{overlap}{right_flank}");
        let r2_sequence = reverse_complement(&oriented_r2_sequence);
        let r1_quality = "I".repeat(r1_sequence.len());
        let r2_quality = "I".repeat(r2_sequence.len());
        let pair = PairInput::new(
            sequence_read(&r1_sequence, &r1_quality),
            sequence_read(&r2_sequence, &r2_quality),
        );
        let mut assembler = exact_match_assembler();

        let search = assembler
            .on_pair(&pair)
            .expect("generated reads should form a valid pair")
            .find_overlap()
            .expect("exact-match overlap search should not error");

        let OverlapOutcome::Found(found) = search else {
            panic!("constructed paired-end insert should produce an overlap");
        };

        let overlap = found.overlap();
        prop_assert_eq!(overlap.len(), overlap_len);
        prop_assert_eq!(overlap.forward_start_offset(), left_flank_len);
        prop_assert_eq!(overlap.reverse_start_offset(), 0);
        prop_assert_eq!(overlap.forward_sequence(), overlap.reverse_sequence());
    }

    #[test]
    fn proptest_incompatible_homopolymer_mates_have_no_exact_overlap(
        r1_len in 20usize..=128,
        r2_len in 20usize..=128,
    ) {
        let r1_sequence = "A".repeat(r1_len);
        let r2_sequence = "C".repeat(r2_len);
        let r1_quality = "I".repeat(r1_sequence.len());
        let r2_quality = "I".repeat(r2_sequence.len());
        let pair = PairInput::new(
            sequence_read(&r1_sequence, &r1_quality),
            sequence_read(&r2_sequence, &r2_quality),
        );
        let mut assembler = exact_match_assembler();

        let search = assembler
            .on_pair(&pair)
            .expect("generated reads should form a valid pair")
            .find_overlap()
            .expect("exact-match overlap search should not error");

        prop_assert!(search.is_no_overlap());
    }
}

fn exact_match_assembler() -> Assembler {
    Assembler::builder()
        .with_overlap_params(OverlapParams::new(0, 20, 0.0, 20))
        .build()
        .expect("exact-match assembler configuration should be valid")
}

fn sequence_read<'a>(sequence: &'a str, quality: &'a str) -> SequenceRead<'a> {
    SequenceRead::try_new("read", sequence, quality)
        .expect("generated read should have matching sequence and quality lengths")
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}
