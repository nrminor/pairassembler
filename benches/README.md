# Comparative real-data benchmarks

This directory contains release-oriented benchmarks that compare `pairasm` with external paired-end merging tools on real ENA data. These are intentionally separate from the Cargo/Criterion benchmarks in `pairassembler/benches/`.

The benchmark harness is a small Cargo package named `pairasm-benches`. It assumes the competitor tools are installed by the benchmark environment rather than by this repository.

Required tools:

- `hyperfine`
- `curl`
- `fastp`
- `bbmerge.sh`
- `vsearch`
- `pairasm`, usually built with `cargo build --release`

The default datasets are listed in `benches/config/datasets.tsv`. Downloads are cached under `benches/data/`, and benchmark runs are written under `benches/runs/`; both locations are ignored by version control.

Run the workflow through `just`:

```sh
just bench-real-check
just bench-real-fetch
READ_PAIRS=100000 just bench-real-prepare
READ_PAIRS=100000 REPLICATES=3 THREADS=8 just bench-real-run
just bench-real-summary
```

You can also call the harness directly:

```sh
cargo run -p pairasm-benches -- check
cargo run -p pairasm-benches -- fetch
cargo run -p pairasm-benches -- prepare --read-pairs 100000
cargo run -p pairasm-benches -- run --read-pairs 100000 --replicates 3 --threads 8
cargo run -p pairasm-benches -- summarize --latest
```

Tool paths can be exported or copied into `benches/config/tools.env` from `tools.env.example`. Benchmark defaults can similarly be copied into `benches/config/benchmark.env` from `benchmark.env.example`.

These comparisons are not perfectly apples-to-apples yet. The tools differ in scoring, filtering, output semantics, and compression behavior. Treat the first runs as decision-support data: runtime, merged counts, and logs should guide which settings and datasets need deeper follow-up before release claims are made.
