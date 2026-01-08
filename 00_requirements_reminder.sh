#!/bin/bash

# Simple reminder of required tools and libraries for PMET/Shiny stack.
# This script does NOT install anything; it only lists what you need and optionally
# tells you whether each binary is currently on PATH.

set -euo pipefail

# Try to load shared color helpers; fall back to plain echo if unavailable.
script_dir=$(cd -- "$(dirname "$0")" && pwd)
if [ -f "$script_dir/scripts/lib/print_colors.sh" ]; then
    # shellcheck source=/dev/null
    source "$script_dir/scripts/lib/print_colors.sh"
else
    print_green() { printf "%s\n" "$1"; }
    print_orange() { printf "%s\n" "$1"; }
    print_fluorescent_yellow() { printf "%s\n" "$1"; }
    print_red() { printf "%s\n" "$1"; }
fi

check_bin() {
    local name="$1"
    if command -v "$name" >/dev/null 2>&1; then
        print_green "[FOUND] $name"
    else
        print_red "[MISSING] $name"
    fi
}

is_build_empty() {
    local build_dir="$script_dir/build"
    if [ ! -d "$build_dir" ]; then
        return 0
    fi
    if find "$build_dir" -mindepth 1 -print -quit >/dev/null 2>&1; then
        return 1
    fi
    return 0
}

setup_pmet_upstream() {
    local external_dir="$script_dir/external"
    local upstream_dir="$external_dir/pmet_project"
    local upstream_repo="https://github.com/duocang/PMET_project"

    if ! command -v git >/dev/null 2>&1; then
        print_red "git is required to fetch upstream PMET sources"
        exit 1
    fi

    mkdir -p "$external_dir"

    if [ ! -d "$upstream_dir/.git" ]; then
        print_orange "[SETUP] Adding PMET_project as submodule at external/pmet_project"
        git submodule add "$upstream_repo" "$upstream_dir"
    else
        print_orange "[SETUP] Updating PMET_project submodule"
        git -C "$upstream_dir" pull --ff-only
    fi

    print_orange "[SETUP] Building upstream PMET binaries"
    (cd "$upstream_dir" && bash scripts/build_all.sh)

    if [ -d "$upstream_dir/build" ]; then
        mkdir -p "$script_dir/build"
        rsync -a --delete "$upstream_dir/build/" "$script_dir/build/"
        print_green "[DONE] Upstream build copied to ./build"
    else
        print_red "[ERROR] Upstream build directory not found after build_all.sh"
        exit 1
    fi
}

maybe_setup_pmet_upstream() {
    if ! is_build_empty; then
        return 0
    fi

    print_orange "[INFO] build/ directory is empty. Run upstream PMET setup? (auto-run in 3s; press n to skip)"
    printf "Run setup now? [Y/n] (auto in 3s): "

    local reply=""
    if read -r -t 3 reply; then
        case "$reply" in
            n|N)
                print_orange "[SKIP] Upstream setup skipped by user"
                return 0
                ;;
        esac
    else
        reply="y"
        printf "\n"
    fi

    if [ -z "$reply" ] || [[ "$reply" =~ ^[Yy]$ ]]; then
        setup_pmet_upstream
    else
        print_orange "[SKIP] Upstream setup skipped"
    fi
}

print_list() {
    local title="$1"; shift
    print_orange "$title"
    printf '%s\n' "$@" | sed 's/^/  - /'
}

# Core binaries
core_bins=(
    "R"
    "Rscript"
    "python3"
    "pip"
    "fasta-get-markov"
    "pmetindex"
    "pmet"
    "pmetParallel"
    "parallel"
    "bedtools"
    "samtools"
)

# R packages (install via scripts/R_utils/install_packages.R)
r_packages=(
    "data.table"
    "tidyverse"
    "ggplot2"
    "hrbrthemes"
    "dplyr"
    "readr"
    "rJava"
)

# Python packages
py_packages=(
    "numpy"
    "pandas"
    "scipy"
    "bio"
    "biopython"
)

print_fluorescent_yellow "PMET / Shiny requirements reminder (no installs performed)"
print_fluorescent_yellow "========================================================"

if [[ ${1-} == "--setup-pmet" ]]; then
    setup_pmet_upstream
    exit 0
fi

# If root build/ is empty, ask (in English) and auto-run setup after 3s unless user declines.
maybe_setup_pmet_upstream

# Binary presence check
print_orange "\nBinary availability (PATH)"
for bin in "${core_bins[@]}"; do
    check_bin "$bin"
done

# Summaries
print_list "\nR packages to have installed" "${r_packages[@]}"
print_list "\nPython packages to have installed" "${py_packages[@]}"

print_orange "\nNotes:"
print_orange "- Install R packages with: Rscript scripts/R_utils/install_packages.R"
print_orange "- Install Python packages with: pip install <pkg>"
print_orange "- MEME/FIMO and PMET binaries are built via scripts/00_binary_compile.sh (interactive)."
print_orange "- Use the built artifacts: copy or symlink external/pmet_project/build to project root if desired"
print_orange "- This script is informational only; it does not change your system."
