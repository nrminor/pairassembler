# pairassembler

[![CI](https://github.com/nrminor/pairassembler/actions/workflows/ci.yml/badge.svg)](https://github.com/nrminor/pairassembler/actions/workflows/ci.yml)
[![Latest release](https://img.shields.io/github/v/release/nrminor/pairassembler?label=release)](https://github.com/nrminor/pairassembler/releases/latest)

The fastest and (laptop-)friendliest paired read merger in the west.

## Overview

`pairassembler` aims to bring overlap-based assembly of paired sequencing reads into the Rust bioinformatics ecosystem, both as a library and as a command line app.

Its main raison d'etre is `libpairassembly`, a Rust library crate for finding no-gap overlaps, validating whether those overlaps are informative, merging mates into a consensus read, and applying overlap-based quality correction. Each of these steps takes inspiration from predecessors including fastp, USEARCH/VSEARCH, and BBMerge, unifying them all into a composable API that tracks each pair's state in the type system.

The `pairassembler` crate, which produces a binary called `pairasm`, dogfoods this library and exposes many of its settings in the command line. It can be thought of as an alternative to BBMerge, fastp, or VSEARCH but solely for pair merging.

## Library quick start

```rust
use libpairassembly::prelude::*;

fn main() -> libpairassembly::Result<()> {
    let forward = SequenceRead::try_new(
        "read-1",
        "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    )?;
    let reverse = SequenceRead::try_new(
        "read-1",
        "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    )?;

    let assembler = Assembler::builder().build()?;
    let pair = PairInput::new(forward, reverse);

    let merged = assembler
        .on_pair(&pair)?
        .find_overlap()?
        .and_then_found(|overlap| overlap.validate()?.merge()?.correct()?.into_owned_read())?;

    if let Some(merged) = merged {
        assert_eq!(merged.id(), "read-1");
        assert_eq!(merged.sequence_bytes().len(), merged.quality_bytes().len());
    }

    Ok(())
}
```

## Installation

A simple curl install is the intended installation path for the command-line tool:

```sh
curl -fsSL https://raw.githubusercontent.com/nrminor/pairassembler/main/INSTALL.sh | bash
```

The installer downloads a pre-built `pairasm` binary for your platform when a release asset is available. If it cannot find a matching binary, it falls back to building from source with Cargo. When a conda, mamba, or pixi-style environment is active, the installer places `pairasm` in that environment's `bin` directory; otherwise it installs to `$HOME/.local/bin`.

## Command line quick start

The CLI is intentionally centered on paired-read merging:

```sh
pairasm \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o merged.fastq.gz \
  --unmerged-out unmerged.fastq.gz
```

Overlap search and validation can be tuned from the command line:

```sh
pairasm \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --min-overlap 30 \
  --overlap-diff-max 2 \
  --diff-percent-max 0.2 \
  --min-complexity-score 39
```

By default, merged reads are quality-corrected using the overlapping evidence from both mates. Use `--no-correct` when you want the merged sequence but do not want overlap-based quality correction.

## What pairassembler does

`pairassembler` treats paired-read assembly as a small, explicit pipeline:

```text
R1/R2 records
    │
    ▼
find ungapped overlap ── no acceptable overlap ──► unmerged branch
    │
    ▼
validate overlap informativeness
    │
    ▼
merge mates into a consensus
    │
    ▼
correct qualities from overlap evidence
```

The main operations are:

- No-gap overlap search between the largest possible overlap windows between mates.
- Validation that separates merely detected overlaps from overlaps informative enough to trust.
- Deterministic consensus merging with configurable tie policies in the library API.
- Overlap-aware quality correction, including a library mode that updates qualities without changing mate bases.
- Normal no-overlap handling: a pair with no acceptable overlap is an expected outcome and not treated as a failure.

## Using the library in another Rust tool

The library does not require applications to adopt a specific FASTQ parser. If your parser record can expose an ID, sequence, and FASTQ ASCII quality string, implement `SeqRecordView` and pass those records through `PairInput`.

The fluent API is useful when you want to make the no-overlap branch explicit or choose a non-default operation order:

```rust
use libpairassembly::prelude::*;

fn merge_pair<R: SeqRecordView>(
    pair: &PairInput<R>,
) -> libpairassembly::Result<Option<OwnedSequenceRead>> {
    let assembler = Assembler::builder()
        .with_overlap_params(OverlapParams::default().with_min_overlap(30))
        .with_validator(OverlapValidator::default().with_min_complexity_score(39))
        .build()?;

    assembler
        .on_pair(pair)?
        .find_overlap()?
        .and_then_found(|overlap| overlap.validate()?.merge()?.correct()?.into_owned_read())
}
```

This shape is deliberately a little more verbose than a one-shot helper. The extra ceremony makes the two runtime branches visible:

```text
find_overlap()? ──► Found(overlap context) ──► validate()? ──► merge()? ──► correct()?
              └──► NoOverlap(context)      ──► Ok(None)
```
