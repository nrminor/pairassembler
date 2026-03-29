# pairassembler project justfile
# All repeating commands should be recipes here.

default:
    @just --list

choose:
    @just --choose

# === Development Workflow ===

check: fmt-check lint test doc-check
    @echo "All checks passed"

check-all: fmt-check lint-strict test-all doc-check
    @echo "All checks passed on full codebase"

# === Formatting ===

fmt-check:
    cargo fmt --all -- --check

fmt:
    cargo fmt --all

# === Linting ===

lint:
    cargo clippy --all-targets --all-features -- -D clippy::correctness -D clippy::unwrap_used

lint-all:
    @just lint

lint-strict:
    cargo clippy --all-targets --all-features -- -D warnings

# === Testing ===

test:
    cargo nextest run --all-features --no-tests=pass

test-all:
    cargo nextest run --all-features --run-ignored all --no-tests=pass

test-verbose:
    cargo nextest run --all-features --no-capture --no-tests=pass

# === Building ===

build:
    cargo build

build-release:
    cargo build --release

check-compile:
    cargo check --all-targets --all-features

# === Documentation ===

doc-check:
    cargo doc --no-deps --document-private-items

doc:
    cargo doc --no-deps --open

# === Project Setup ===

setup: clone-refs
    @echo ""
    @echo "Project setup complete"
    @echo "Reference repos: .agents/repos/"

# === Reference Repositories ===

[arg("force", long="force", value="1")]
clone-refs force="":
    @echo "Cloning reference repositories into .agents/repos/..."
    @{{ if force != "" { "rm -rf .agents/repos" } else { "true" } }}
    @mkdir -p .agents/repos
    @echo "Cloning vsearch..."
    git clone --depth 1 https://github.com/torognes/vsearch.git .agents/repos/vsearch || echo "vsearch already exists, skipping"
    @echo "Cloning fastp..."
    git clone --depth 1 https://github.com/OpenGene/fastp.git .agents/repos/fastp || echo "fastp already exists, skipping"
    @echo "Cloning bbmap..."
    git clone --depth 1 https://github.com/BioInfoTools/BBMap.git .agents/repos/bbmap || echo "bbmap already exists, skipping"
    @echo "Cloning FLASH2..."
    git clone --depth 1 https://github.com/dstreett/FLASH2.git .agents/repos/flash2 || echo "flash2 already exists, skipping"
    @echo "Cloning SeqPrep..."
    git clone --depth 1 https://github.com/jstjohn/SeqPrep.git .agents/repos/seqprep || echo "seqprep already exists, skipping"
    @echo "Cloning PEAR..."
    git clone --depth 1 https://github.com/tseemann/PEAR.git .agents/repos/pear || echo "pear already exists, skipping"
    @echo "Cloning NGmerge..."
    git clone --depth 1 https://github.com/jsh58/NGmerge.git .agents/repos/ngmerge || echo "ngmerge already exists, skipping"
    @echo "Reference repositories cloned to .agents/repos/"

clean-refs:
    @echo "Removing reference repositories..."
    rm -rf .agents/repos
    @echo "Reference repositories removed"

# === Utility ===

clean:
    cargo clean

update:
    cargo update

sloc:
    @tokei --types=Rust --compact
