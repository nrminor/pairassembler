use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use libpairassembly::{Assembler, OverlapSearch, PairInput, SequenceRead};

#[derive(Clone, Copy, Debug)]
enum PairKind {
    Mergeable,
    NoOverlap,
    ValidationRejected,
}

struct BenchPair {
    id: String,
    r1_sequence: String,
    r2_sequence: String,
    quality: String,
}

fn compute_mixed_pairs(c: &mut Criterion) {
    let assembler = Assembler::builder().build().unwrap_or_else(|error| {
        panic!("failed to build default assembler for benchmark: {error}");
    });
    let pairs = build_compute_pairs();

    c.bench_function("compute_mixed_pairs_300", |b| {
        b.iter(|| count_merged_pairs(&assembler, &pairs));
    });
}

fn build_compute_pairs() -> Vec<BenchPair> {
    let mut pairs = Vec::new();
    pairs.extend(many_pairs(100, PairKind::Mergeable));
    pairs.extend(many_pairs(100, PairKind::NoOverlap));
    pairs.extend(many_pairs(100, PairKind::ValidationRejected));
    pairs
}

fn count_merged_pairs(assembler: &Assembler, pairs: &[BenchPair]) -> usize {
    let mut merged = 0_usize;

    for pair in pairs {
        let read1 = read_from_parts(&pair.id, &pair.r1_sequence, &pair.quality);
        let read2 = read_from_parts(&pair.id, &pair.r2_sequence, &pair.quality);
        let input = PairInput::new(read1, read2);
        if assemble_pair(assembler, &input) {
            merged = merged.saturating_add(1);
        }
    }

    black_box(merged)
}

fn assemble_pair(assembler: &Assembler, input: &PairInput<SequenceRead<'_>>) -> bool {
    let ready = assembler
        .on_pair(input)
        .unwrap_or_else(|error| panic!("unexpected pair setup error in benchmark: {error}"));
    let overlap = match ready
        .find_overlap()
        .unwrap_or_else(|error| panic!("unexpected overlap error in benchmark: {error}"))
    {
        OverlapSearch::Found(overlap) => overlap,
        OverlapSearch::NoOverlap(_) => return false,
    };
    let validated = match overlap.validate() {
        Ok(validated) => validated,
        Err(_validation_rejected) => return false,
    };
    let merged = validated
        .merge()
        .unwrap_or_else(|error| panic!("unexpected merge error in benchmark: {error}"));
    let corrected = merged
        .correct()
        .unwrap_or_else(|error| panic!("unexpected correction error in benchmark: {error}"));
    let _read = corrected.into_owned_read().unwrap_or_else(|error| {
        panic!("unexpected owned-read conversion error in benchmark: {error}")
    });
    true
}

fn read_from_parts<'a>(id: &'a str, sequence: &'a str, quality: &'a str) -> SequenceRead<'a> {
    SequenceRead::try_new(id, sequence, quality).unwrap_or_else(|error| {
        panic!("invalid benchmark read fixture: {error}");
    })
}

fn many_pairs(count: usize, kind: PairKind) -> Vec<BenchPair> {
    (0..count)
        .map(|index| pair(&format!("{kind:?}-{index:08}"), kind))
        .collect()
}

fn pair(id: &str, kind: PairKind) -> BenchPair {
    match kind {
        PairKind::Mergeable => BenchPair::new(
            id,
            "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
            "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
        ),
        PairKind::NoOverlap => BenchPair::homopolymer(id, 'A', 'C', 56),
        PairKind::ValidationRejected => BenchPair::homopolymer(id, 'A', 'T', 56),
    }
}

impl BenchPair {
    fn new(id: &str, r1_sequence: &str, r2_sequence: &str) -> Self {
        Self {
            id: id.to_owned(),
            r1_sequence: r1_sequence.to_owned(),
            r2_sequence: r2_sequence.to_owned(),
            quality: "I".repeat(r1_sequence.len()),
        }
    }

    fn homopolymer(id: &str, r1_base: char, r2_base: char, len: usize) -> Self {
        Self {
            id: id.to_owned(),
            r1_sequence: r1_base.to_string().repeat(len),
            r2_sequence: r2_base.to_string().repeat(len),
            quality: "I".repeat(len),
        }
    }
}

criterion_group!(benches, compute_mixed_pairs);
criterion_main!(benches);
