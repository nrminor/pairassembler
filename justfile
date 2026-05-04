# pairassembler project justfile.
# All repeating commands should be recipes here.

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

# Run quick benchmark smoke checks.
bench-smoke:
    cargo bench --bench compute -- --test
    PAIRASM_E2E_PAIRS=1000 cargo bench --bench e2e_pipeline -- --test

# Run the compute-only Criterion benchmark.
bench-compute:
    cargo bench --bench compute

# Run the synthetic e2e Criterion benchmark.
bench-e2e:
    PAIRASM_E2E_PAIRS=10000 cargo bench --bench e2e_pipeline

# Check external tools for real-data comparative benchmarks.
bench-real-check:
    cargo run -p pairasm-benches -- check

# Fetch configured ENA datasets for real-data comparative benchmarks.
bench-real-fetch:
    cargo run -p pairasm-benches -- fetch --config benches/config/datasets.tsv

# Prepare deterministic first-N-pair subsets for real-data benchmarks.
bench-real-prepare:
    cargo run -p pairasm-benches -- prepare --config benches/config/datasets.tsv --read-pairs ${READ_PAIRS:-100000}

# Run pairasm and competitor real-data benchmarks through hyperfine.
bench-real-run: build-release
    PAIRASM_BIN=${PAIRASM_BIN:-target/release/pairasm} cargo run -p pairasm-benches -- run --config benches/config/datasets.tsv --read-pairs ${READ_PAIRS:-100000} --replicates ${REPLICATES:-3} --threads ${THREADS:-8}

# Summarize the latest real-data benchmark run.
bench-real-summary:
    cargo run -p pairasm-benches -- summarize --latest

alias bs := bench-smoke
alias bc := bench-compute
alias be := bench-e2e
alias brc := bench-real-check
alias brf := bench-real-fetch
alias brp := bench-real-prepare
alias brr := bench-real-run
alias brs := bench-real-summary

# === Building ===

# Build the workspace in debug mode.
build:
    cargo build

# Build the workspace in release mode.
build-release:
    cargo build --release

# Typecheck all targets and features.
check-compile:
    cargo check --all-targets --all-features

alias b := build
alias r := build-release
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
