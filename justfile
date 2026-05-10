# pairassembler project justfile.
# All repeating commands should be recipes here.

export DUCKDB_DOWNLOAD_LIB := env_var_or_default("DUCKDB_DOWNLOAD_LIB", "1")

# Show the available recipes by default.
default:
    @just --list

# Interactively choose a recipe to run.
choose:
    @just --choose

# List all available recipes.
list:
    @just --list

alias l := list
alias ls := list
alias ch := choose

# === Development Workflow ===

# Run the regular local check suite.
check: fmt-check lint test doc-check
    @echo "All checks passed"

# Run the strict full-codebase check suite.
check-all: fmt-check lint-strict test-all doc-check
    @echo "All checks passed on full codebase"

alias c := check
alias ca := check-all

# Prepare a commit by running the local check suite and showing jj status.
prepare-commit: check
    @echo ""
    @echo "Ready to describe or split the current jj commit."
    @jj status

# Prepare for push by running the strict full-codebase check suite.
prepare-push: check-all
    @echo ""
    @echo "Ready to push when you are."

alias pc := prepare-commit
alias pp := prepare-push

# === Formatting ===

# Check formatting without modifying files.
fmt-check:
    cargo fmt --all -- --check

# Apply formatting fixes.
fmt:
    cargo fmt --all

alias fc := fmt-check
alias f := fmt

# === Linting ===

# Run the fast local Clippy policy.
lint:
    cargo clippy --all-targets --all-features -- -D clippy::correctness -D clippy::unwrap_used

# Compatibility alias for the standard lint recipe name.
lint-all:
    @just lint

# Run strict Clippy with all warnings denied.
lint-strict:
    cargo clippy --all-targets --all-features -- -D warnings

alias cl := lint
alias cla := lint-strict

# === Testing ===

# Run tests with nextest.
test:
    cargo nextest run --all-features --no-tests=pass

# Run all tests, including ignored tests.
test-all:
    cargo nextest run --all-features --run-ignored all --no-tests=pass

# Run tests with captured output disabled.
test-verbose:
    cargo nextest run --all-features --no-capture --no-tests=pass

alias tr := test
alias ta := test-all
alias tv := test-verbose

# === Benchmarks ===

# Fast local pairasm performance sanity check.
benchmark-smoke: _benchmark-pairasm-smoke

# Run the standard real-data tool comparison and print a report.
benchmark: _benchmark-default

# Run the standard real-data tool comparison and print a report.
_benchmark-default: _benchmark-check-tools _benchmark-fetch-ena _benchmark-prepare-subsets _benchmark-run-default _benchmark-report-read-id-overlap

# Run the tuned/comparability tool comparison and print a report.
benchmark-tuned: _benchmark-check-tools _benchmark-fetch-ena _benchmark-prepare-subsets _benchmark-run-tuned _benchmark-report-tuned-read-id-overlap

# Reprint the latest default-user benchmark comparison report.
benchmark-report: _benchmark-report-read-id-overlap

# Verify Criterion benchmark targets run quickly; this is not a measurement.
_benchmark-pairasm-smoke:
    @just _benchmark-phase "smoke" "Checking benchmark targets compile and run quickly"
    @cargo bench --bench in_memory_merge -- --test
    @PAIRASM_FASTQ_PAIRS=1000 cargo bench --bench fastq_merge -- --test

# Measure pairasm's in-memory merge path with Criterion.
_benchmark-pairasm-in-memory:
    cargo bench --bench in_memory_merge

# Measure pairasm's synthetic FASTQ-oriented path with Criterion.
_benchmark-pairasm-fastq:
    PAIRASM_FASTQ_PAIRS=10000 cargo bench --bench fastq_merge

# Check external tools needed for pairasm-vs-tool comparison runs.
_benchmark-check-tools:
    @just _benchmark-phase "1/5" "Checking external benchmark tools"
    @cargo run --quiet -p pairasm-benches -- check

# Fetch configured ENA FASTQ inputs for pairasm-vs-tool comparisons.
_benchmark-fetch-ena:
    @just _benchmark-phase "2/5" "Fetching configured ENA FASTQs"
    @cargo run --quiet -p pairasm-benches -- fetch --config benches/config/datasets.tsv

# Prepare deterministic first-N-pair subsets from fetched ENA inputs.
_benchmark-prepare-subsets:
    @just _benchmark-phase "3/5" "Preparing deterministic FASTQ subsets"
    @cargo run --quiet -p pairasm-benches -- prepare --config benches/config/datasets.tsv --read-pairs ${READ_PAIRS:-100000}

# Run the default-user pairasm-vs-tool comparison through hyperfine.
_benchmark-run-default:
    @just _benchmark-phase "4/5" "Running default-user tool comparison"
    @cargo run --quiet -p pairasm-benches -- run --config benches/config/datasets.tsv --read-pairs ${READ_PAIRS:-100000} --replicates ${REPLICATES:-3} --threads ${THREADS:-8} --mode default-user --db ${BENCHMARK_DB:-benches/benchmarks.duckdb}

# Run the tuned/comparability pairasm-vs-tool comparison through hyperfine.
_benchmark-run-tuned:
    @just _benchmark-phase "4/5" "Running tuned/comparability tool comparison"
    @cargo run --quiet -p pairasm-benches -- run --config benches/config/datasets.tsv --read-pairs ${READ_PAIRS:-100000} --replicates ${REPLICATES:-3} --threads ${THREADS:-8} --mode tuned-comparability --db ${BENCHMARK_DB:-benches/benchmarks.duckdb}

# Report merged read-ID set overlap for the latest default-user benchmark run in DuckDB.
_benchmark-report-read-id-overlap:
    @just _benchmark-phase "5/5" "Reporting merged read-ID set overlap"
    @cargo run --quiet -p pairasm-benches -- report read-id-overlap --db ${BENCHMARK_DB:-benches/benchmarks.duckdb}

# Report merged read-ID set overlap for the latest tuned benchmark run in DuckDB.
_benchmark-report-tuned-read-id-overlap:
    @just _benchmark-phase "5/5" "Reporting tuned merged read-ID set overlap"
    @cargo run --quiet -p pairasm-benches -- report read-id-overlap --mode tuned-comparability --db ${BENCHMARK_DB:-benches/benchmarks.duckdb}

# Print a compact visual separator for benchmark phases.
_benchmark-phase step title:
    @cargo run --quiet -p pairasm-benches -- workflow-phase "{{step}}" "{{title}}"

alias bm := benchmark

# === Building ===

# Build the workspace in debug mode.
build:
    cargo build

# Build the workspace in release mode.
build-release:
    cargo build --release

# Install the `pairasm` binary globally on the user's system
install:
    cargo install --path ./pairassembler/

# Typecheck all targets and features.
check-compile:
    cargo check --all-targets --all-features

# Generate the experimental C ABI header with cbindgen.
c-header:
    mkdir -p c-pairassembler/include
    cbindgen --config cbindgen.toml --crate libpairassembly --output c-pairassembler/include/libpairassembly.h

alias b := build
alias r := build-release
alias i := install
alias cr := check-compile

# === Documentation ===

# Check that documentation builds without errors.
doc-check:
    cargo doc --no-deps --document-private-items

# Generate and open documentation.
doc:
    cargo doc --no-deps --open

alias dc := doc-check
alias d := doc

# === Project Setup ===

# Prepare reference repositories for manual comparison work.
setup: clone-refs
    @echo ""
    @echo "Project setup complete"
    @echo "Reference repos: .agents/repos/"

# === Reference Repositories ===

# Clone external reference implementations into .agents/repos.
[arg("force", long="force", value="1")]
clone-refs force="":
    @echo "Cloning reference repositories into .agents/repos/..."
    @{{ if force != "" { "rm -rf .agents/repos" } else { "true" } }}
    @mkdir -p .agents/repos
    @mkdir -p .agents/archive
    @echo "Cloning vsearch..."
    git clone --depth 1 https://github.com/torognes/vsearch.git .agents/repos/vsearch || echo "vsearch already exists, skipping"
    @echo "Cloning fastp..."
    git clone --depth 1 https://github.com/OpenGene/fastp.git .agents/repos/fastp || echo "fastp already exists, skipping"
    @echo "Downloading bbmap from SourceForge..."
    @if [ -d .agents/repos/bbmap ]; then echo "bbmap already exists, skipping"; else curl -L "https://sourceforge.net/projects/bbmap/files/latest/download" -o .agents/archive/bbmap-latest.tar.gz && tar -xzf .agents/archive/bbmap-latest.tar.gz -C .agents/repos; fi
    @echo "Cloning FLASH2..."
    git clone --depth 1 https://github.com/dstreett/FLASH2.git .agents/repos/flash2 || echo "flash2 already exists, skipping"
    @echo "Cloning SeqPrep..."
    git clone --depth 1 https://github.com/jstjohn/SeqPrep.git .agents/repos/seqprep || echo "seqprep already exists, skipping"
    @echo "Cloning PEAR..."
    git clone --depth 1 https://github.com/tseemann/PEAR.git .agents/repos/pear || echo "pear already exists, skipping"
    @echo "Cloning NGmerge..."
    git clone --depth 1 https://github.com/jsh58/NGmerge.git .agents/repos/ngmerge || echo "ngmerge already exists, skipping"
    @echo "Reference repositories cloned to .agents/repos/"

# Remove downloaded reference repositories.
clean-refs:
    @echo "Removing reference repositories..."
    rm -rf .agents/repos
    @echo "Reference repositories removed"

# === Version Control ===

# Show current jj status.
status:
    jj status

# Show current jj log.
log:
    jj log

alias st := status
alias lg := log

# === Utility ===

# Clean build artifacts.
clean:
    cargo clean

# Update dependencies.
update:
    cargo update

# Count source lines.
sloc:
    @tokei --types=Rust --compact
