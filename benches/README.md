# Comparative real-data benchmarks

This directory contains release-oriented benchmarks that compare `pairasm` with external paired-end merging tools on real ENA data. These are intentionally separate from the Cargo/Criterion benchmarks in `pairassembler/benches/`.

The benchmark harness is a small Cargo package named `pairasm-benches`. It assumes the competitor tools are installed by the benchmark environment rather than by this repository.

Required tools:

- `hyperfine`
- `curl`
- `fastp`
- `bbmerge` or legacy `bbmerge.sh`
- `vsearch`
- `pairasm`, usually built with `cargo build --release`

The default datasets are listed in `benches/config/datasets.tsv`. Downloads are cached under `benches/data/`, and benchmark evidence is written under `benches/runs/`; both locations are ignored by version control.

Run the workflow through `just`:

```sh
just bench
READ_PAIRS=100000 REPLICATES=3 THREADS=8 just benchmark
```

`bench` is a fast local sanity check for pairasm's Criterion benchmarks. `benchmark` is the standard end-to-end real-data comparison: it builds the release binary, checks external tools, fetches missing ENA inputs, prepares deterministic subsets, runs each merge tool, validates the outputs, and prints a comparison report.

The standard comparison uses the `default-user` mode, which is intended to model a hurried user running each tool directly from paired R1/R2 FASTQs with minimal extra thought. A tuned/comparability mode is available when investigating how tools behave under closer merge policies:

```sh
READ_PAIRS=100000 REPLICATES=3 THREADS=8 just benchmark-tuned
```

To reprint the latest comparison report without running a benchmark, use `just benchmark-report`.

You can also call the harness directly:

```sh
cargo run -p pairasm-benches -- check
cargo run -p pairasm-benches -- fetch
cargo run -p pairasm-benches -- prepare --read-pairs 100000
cargo run -p pairasm-benches -- run --read-pairs 100000 --replicates 3 --threads 8 --mode default-user
cargo run -p pairasm-benches -- report agreement
```

Tool paths can be exported or copied into `benches/config/tools.env` from `tools.env.example`. Benchmark defaults can similarly be copied into `benches/config/benchmark.env` from `benchmark.env.example`.

Structured benchmark results are stored in `benches/benchmarks.duckdb` so reports can be regenerated without rerunning the tools. Set `BENCHMARK_DB=/path/to/benchmarks.duckdb` if you need to read or write a non-default results store.

Preparation writes paired R1/R2 subsets. All comparison tools run directly from those paired FASTQs in both benchmark modes.

The `default-user` mode leaves VSEARCH's merge policy at its CLI defaults, apart from required input/output paths and the requested thread count. That conservatism is part of what the benchmark is meant to reveal. The `tuned-comparability` mode instead allows staggered merges, sets the minimum overlap to 30 bases, allows at most 5 differences, and caps the overlap difference percentage at 20%. Those settings are useful for investigation but should not be confused with the minimal-thought default benchmark.

The run step validates tool accounting where the external tools report it. In particular, BBMerge and VSEARCH must report the expected number of input pairs, their merged output count must match the output FASTQ, and BBMerge stderr must not contain Java exceptions. This is intentionally strict: a benchmark artifact that only partially represents a tool run degrades trust more than a failed benchmark does.

It should be emphasized that these comparisons are knowingly not perfectly apples-to-apples. The tools differ in scoring, filtering, output semantics, and compression behavior. In particular, fastp's default merge output reflects both overlap detection and fastp's general read filtering, while `pairasm` is focused on overlap-based consensus assembly and does not apply the same pre-merge read-level filters. This makes default CLI comparisons useful for understanding normal tool behavior. They should not, however, be taken as a pure comparison of overlap discovery.

A `pairasm`-only merged read may therefore represent an intended salvage case: a pair whose individual reads look low quality, but whose overlapping evidence supports a corrected consensus. Treat benchmark results as decision-support data. Runtime, merged counts, output membership, and logs should guide which settings and datasets need deeper follow-up before release claims are made.
