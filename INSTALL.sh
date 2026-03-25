#!/usr/bin/env bash
set -euo pipefail

# pairassembler curl installer

REPO="nrminor/pairassembler"
BINARY_NAME="pairasm"
VERSION="0.1.0"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
step() { echo -e "${BLUE}[STEP]${NC} $1"; }

in_conda_env() {
	[[ -n "${CONDA_PREFIX:-}" ]]
}

determine_install_dir() {
	if in_conda_env; then
		echo "${CONDA_PREFIX}/bin"
	else
		echo "${HOME}/.local/bin"
	fi
}

get_cargo_install_root() {
	local install_dir="$1"
	echo "${install_dir%/bin}"
}

get_environment_name() {
	if [[ -n "${PIXI_PROJECT_ROOT:-}" ]]; then
		echo "pixi"
	elif [[ -n "${CONDA_PREFIX:-}" ]]; then
		if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
			echo "conda (${CONDA_DEFAULT_ENV})"
		else
			echo "conda"
		fi
	else
		echo "none"
	fi
}

detect_platform() {
	local os arch

	case "$(uname -s)" in
	Linux) os="linux" ;;
	Darwin) os="darwin" ;;
	MINGW* | MSYS* | CYGWIN*) os="windows" ;;
	*)
		error "Unsupported OS: $(uname -s)"
		exit 1
		;;
	esac

	case "$(uname -m)" in
	x86_64 | amd64) arch="x86_64" ;;
	aarch64 | arm64) arch="aarch64" ;;
	*)
		error "Unsupported architecture: $(uname -m)"
		exit 1
		;;
	esac

	case "${os}-${arch}" in
	linux-x86_64) echo "x86_64-unknown-linux-musl" ;;
	linux-aarch64) echo "aarch64-unknown-linux-musl" ;;
	darwin-x86_64) echo "x86_64-apple-darwin" ;;
	darwin-aarch64) echo "aarch64-apple-darwin" ;;
	windows-x86_64) echo "x86_64-pc-windows-msvc" ;;
	windows-aarch64) echo "aarch64-pc-windows-msvc" ;;
	*)
		error "Unsupported platform: ${os}-${arch}"
		exit 1
		;;
	esac
}

get_latest_release() {
	local url="https://api.github.com/repos/${REPO}/releases/latest"
	if command -v curl &>/dev/null; then
		curl -fsSL "$url" | grep '"tag_name":' | sed -E 's/.*"([^"]+)".*/\1/'
	elif command -v wget &>/dev/null; then
		wget -qO- "$url" | grep '"tag_name":' | sed -E 's/.*"([^"]+)".*/\1/'
	else
		error "Neither curl nor wget found. Please install one of them."
		exit 1
	fi
}

install_binary() {
	local version="$1"
	local platform="$2"
	local install_dir="$3"

	local archive_ext="tar.gz"
	if [[ "$platform" == *"windows"* ]]; then
		archive_ext="zip"
	fi

	local download_url="https://github.com/${REPO}/releases/download/${version}/${BINARY_NAME}-${platform}.${archive_ext}"
	local temp_dir
	temp_dir=$(mktemp -d)

	step "Downloading ${BINARY_NAME} ${version} for ${platform}..."
	cd "$temp_dir"

	if command -v curl &>/dev/null; then
		curl -fsSL -o "archive.${archive_ext}" "$download_url" || return 1
	else
		wget -q -O "archive.${archive_ext}" "$download_url" || return 1
	fi

	step "Extracting binary..."
	if [[ "$archive_ext" == "zip" ]]; then
		unzip -q "archive.${archive_ext}" || return 1
	else
		tar -xzf "archive.${archive_ext}" || return 1
	fi

	step "Installing to ${install_dir}..."
	mkdir -p "$install_dir"

	if [[ -f "${BINARY_NAME}" ]]; then
		chmod +x "${BINARY_NAME}"
		mv "${BINARY_NAME}" "${install_dir}/"
	elif [[ -f "${BINARY_NAME}.exe" ]]; then
		mv "${BINARY_NAME}.exe" "${install_dir}/"
	else
		error "Binary not found in archive"
		return 1
	fi

	cd - >/dev/null || true
	rm -rf "$temp_dir"

	info "Successfully installed ${BINARY_NAME} to ${install_dir}"
}

fail_no_cargo_in_env() {
	local env_type="$1"

	error "Cargo not found in active environment"
	echo ""

	case "$env_type" in
	pixi*)
		info "Add Rust to your pixi.toml:"
		echo "  pixi add rust"
		;;
	conda*)
		info "Install Rust in your conda/mamba environment:"
		echo "  conda install -c conda-forge rust"
		echo "  # or"
		echo "  mamba install rust"
		;;
	*)
		info "Install Rust globally:"
		echo "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
		;;
	esac

	exit 1
}

build_from_source() {
	local install_dir="$1"
	local install_root
	install_root=$(get_cargo_install_root "$install_dir")

	step "Building from source..."

	if ! command -v cargo &>/dev/null; then
		if in_conda_env; then
			local env_type
			env_type=$(get_environment_name)
			fail_no_cargo_in_env "$env_type"
		else
			warn "Rust not found. Installing Rust toolchain..."
			curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
			export PATH="$HOME/.cargo/bin:$PATH"
		fi
	fi

	info "Installing ${BINARY_NAME} from GitHub repository..."
	if cargo install --git "https://github.com/${REPO}.git" --root "${install_root}" --force pairassembler --bin pairasm; then
		info "Successfully built and installed ${BINARY_NAME}"
	else
		error "Failed to install from source"
		exit 1
	fi
}

show_help() {
	echo -e "${BOLD}pairassembler Installer${NC} - Version ${VERSION}"
	echo ""
	echo "Installs the pairasm CLI either from prebuilt release assets or from source."
	echo ""
	echo "Usage: $0 [--help] [--version]"
}

show_version() {
	cat <<EOF
${BINARY_NAME} installer version ${VERSION}
Repository: https://github.com/${REPO}
License: MIT
EOF
}

parse_args() {
	while [[ $# -gt 0 ]]; do
		case $1 in
		-h | --help)
			show_help
			exit 0
			;;
		-v | --version)
			show_version
			exit 0
			;;
		*)
			error "Unknown option: $1"
			echo "Use --help for usage information"
			exit 1
			;;
		esac
	done
}

verify_installation() {
	local install_dir="$1"

	if command -v "${BINARY_NAME}" &>/dev/null; then
		echo -e "${GREEN}${BOLD}✓ Installation successful!${NC}"
		info "Run '${BINARY_NAME} --help' to get started"
	else
		warn "Binary installed but not found in PATH"
		info "Run directly: ${install_dir}/${BINARY_NAME}"
	fi
}

main() {
	parse_args "$@"

	echo ""
	echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	echo -e "${BOLD}  pairassembler Installer${NC}"
	echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	echo ""

	local install_dir
	install_dir=$(determine_install_dir)
	info "Install directory: ${install_dir}"

	local platform
	platform=$(detect_platform)
	info "Detected platform: ${platform}"

	local version
	version=$(get_latest_release || true)

	if [[ -n "${version:-}" ]]; then
		info "Latest release: ${version}"
		if install_binary "$version" "$platform" "$install_dir"; then
			verify_installation "$install_dir"
			return 0
		fi
		warn "Failed to download pre-built binary"
	else
		warn "Could not determine latest release"
	fi

	warn "Falling back to source build..."
	build_from_source "$install_dir"
	verify_installation "$install_dir"
}

main "$@"
