#!/bin/bash
# ==============================================================================
# PMET Index for Genomic Intervals
# ==============================================================================
# Author  : Charlotte Rich (22.1.18)
# Purpose : Creates PMET index for Paired Motif Enrichment Test
# Requires: scripts/lib/print_colors.sh, scripts/lib/timer.sh
# ==============================================================================

set -e
set -o errexit
set -o pipefail

# ==============================================================================
# INITIALIZATION
# ==============================================================================

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
pmetDir="build"

# Source shared libraries
for lib in print_colors timer; do
    lib_file="$project_root/scripts/lib/$lib.sh"
    if [ ! -f "$lib_file" ]; then
        echo "ERROR: Missing library at $lib_file" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$lib_file"
done

start_time=$SECONDS

# ==============================================================================
# FUNCTIONS
# ==============================================================================

function usage() {
    cat >&2 <<EOF

USAGE: $(basename "$0") [options] <genome> <memefile>

Creates PMET index for Paired Motif Enrichment Test using genomic intervals.

POSITIONAL ARGUMENTS:
    <genome>              Genomic interval FASTA file
    <memefile>            MEME format motif file

REQUIRED OPTIONS:
    -o <directory>        Output directory for results

OPTIONAL ARGUMENTS:
    -r <path>             Scripts root directory         [default: scripts]
    -n <topn>             Top N hits per motif           [default: 5000]
    -k <maxk>             Max K hits per interval        [default: 5]
    -f <threshold>        FIMO p-value threshold         [default: 0.05]
    -t <threads>          Parallel threads               [default: 4]
    -h                    Show this help message

EXAMPLE:
    $(basename "$0") -o results/output -t 8 genome.fa motifs.meme

EOF
}

function error_exit() {
    print_red "\n[ERROR] $1"
    usage
    print_elapsed_time_with_status "$start_time" "error" "$1"
    exit 1
}

# ==============================================================================
# DEFAULT ARGUMENTS
# ==============================================================================

topn=5000
maxk=5
fimothresh=0.05
pmetroot="scripts"
threads=4

outputdir=""
genomefile=""
memefile=""

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

if [ $# -eq 0 ]; then
    error_exit "No arguments supplied"
fi

while getopts ":r:o:k:n:f:t:h" opt; do
    case $opt in
        r) pmetroot="$OPTARG"   ;;
        o) outputdir="$OPTARG"  ;;
        n) topn="$OPTARG"       ;;
        k) maxk="$OPTARG"       ;;
        f) fimothresh="$OPTARG" ;;
        t) threads="$OPTARG"    ;;
        h) usage; exit 0        ;;
        \?) error_exit "Invalid option: -$OPTARG" ;;
        :)  error_exit "Option -$OPTARG requires an argument" ;;
    esac
done

shift $((OPTIND - 1))
genomefile="${1:-}"
memefile="${2:-}"

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

[ -z "$genomefile" ] && error_exit "Missing required argument: <genome> file"
[ -z "$memefile" ]   && error_exit "Missing required argument: <memefile>"
[ -z "$outputdir" ]  && error_exit "Missing required option: -o <output_directory>"
[ ! -f "$genomefile" ] && error_exit "Genome file not found: $genomefile"
[ ! -f "$memefile" ]   && error_exit "Meme file not found: $memefile"

# ==============================================================================
# DISPLAY CONFIGURATION
# ==============================================================================

print_green "\n╔══════════════════════════════════════════════════════════════════╗"
print_green "║                    PMET INDEX CONFIGURATION                      ║"
print_green "╚══════════════════════════════════════════════════════════════════╝"
print_white "  Genome file            : "; print_orange "$genomefile"
print_white "  Meme file              : "; print_orange "$memefile"
print_white "  Output directory       : "; print_orange "$outputdir"
print_green "────────────────────────────────────────────────────────────────────"
print_white "  Top N hits per motif   : "; print_orange "$topn"
print_white "  Max K hits per interval: "; print_orange "$maxk"
print_white "  FIMO threshold         : "; print_orange "$fimothresh"
print_white "  Parallel threads       : "; print_orange "$threads"
print_green "────────────────────────────────────────────────────────────────────\n"

mkdir -p "$outputdir"

# ==============================================================================
# STEP 1: PREPARE SEQUENCES
# ==============================================================================

# Preprocess FASTA - replace ':' with '__COLON__' in sequence names
# FIMO has issues parsing sequence names containing ':'
genomefile_temp="${genomefile%.fa}_temp.fa"
sed 's/^\(>.*\):/\1__COLON__/g' "$genomefile" > "$genomefile_temp"
print_orange "   └─ Created temporary file with sanitized sequence names"



print_green "┌─ Step 1/4: Preparing sequences..."
step_start=$SECONDS

universefile="$outputdir/universe.txt"

if [[ ! -f "$universefile" || ! -f "$outputdir/promoter_lengths.txt" ]]; then
    python3 "$pmetroot/python/deduplicate.py" \
            "$genomefile_temp" \
            "$outputdir/no_duplicates.fa" \
        || error_exit "Deduplication failed"

    python3 "$pmetroot/python/parse_promoter_lengths_from_fasta.py" \
            "$outputdir/no_duplicates.fa" \
            "$outputdir/promoter_lengths.txt" \
        || error_exit "Promoter length parsing failed"

    cut -f1 "$outputdir/promoter_lengths.txt" > "$universefile"
    [ -f "$outputdir/no_duplicates.fa" ] && rm -f "$outputdir/no_duplicates.fa"
fi

print_elapsed_time "$step_start"

# ==============================================================================
# STEP 2: BUILD BACKGROUND MODEL & PROCESS MOTIFS
# ==============================================================================

print_green "┌─ Step 2/4: Building background model and processing motifs..."
step_start=$SECONDS

fasta-get-markov "$genomefile_temp" > "$outputdir/genome.bg" \
    || error_exit "fasta-get-markov failed"

mkdir -p "$outputdir/memefiles"

python3 "$pmetroot/python/parse_memefile.py" \
        "$memefile" \
        "$outputdir/memefiles/" \
    || error_exit "parse_memefile failed"

python3 "$pmetroot/python/calculateICfrommeme_IC_to_csv.py" \
        "$outputdir/memefiles/" \
        "$outputdir/IC.txt" \
    || error_exit "calculateICfrommeme failed"

mkdir -p "$outputdir/fimo" "$outputdir/fimohits"

shopt -s nullglob
nummotifs=$(grep -c '^MOTIF' "$memefile")
print_orange "   └─ Found $nummotifs motifs"

print_elapsed_time "$step_start"

# # ==============================================================================
# # STEP 3: RUN FIMO IN PARALLEL
# # ==============================================================================

# print_green "┌─ Step 3/4: Running FIMO ($threads parallel threads)..."
# step_start=$SECONDS

# n=0
# for meme_file in "$outputdir/memefiles"/*.txt; do
#     ((n++))
#     fimofile=$(basename "$meme_file")
#     fimo --no-qvalue --text --thresh "$fimothresh" --verbosity 1 \
#          --bgfile "$outputdir/genome.bg" "$meme_file" "$genomefile_temp" \
#          > "$outputdir/fimo/$fimofile" &
#     [ $((n % threads)) -eq 0 ] && wait
# done
# wait

# # Post-process FIMO output - restore ':' in sequence names
# for fimofile in "$outputdir/fimo"/*.txt; do
#     sed -i '' 's/__COLON__/:/g' "$fimofile"
# done
# rm -f "$genomefile_temp"
# print_orange "   └─ Restored sequence names and cleaned up temporary file"

# print_elapsed_time "$step_start"

# # ==============================================================================
# # STEP 4: INDEX AND FILTER MOTIF HITS
# # ==============================================================================

# print_green "┌─ Step 4/4: Indexing and filtering motif hits..."
# step_start=$SECONDS

# "$pmetDir/index_cpp" \
#     -f "$outputdir/fimo" \
#     -k "$maxk" \
#     -n "$topn" \
#     -p "$outputdir/promoter_lengths.txt" \
#     -o "$outputdir" \
#     || error_exit "pmetindex failed"

# print_elapsed_time "$step_start"

n=0
for meme_file in "$outputdir/memefiles"/*.txt; do
    ((n++))
    fimofile=$(basename "$meme_file")
    "$pmetDir/index_fimo_fused"             \
        --no-qvalue                         \
        --text                              \
        --thresh "$fimothresh"              \
        --verbosity 1                       \
        --bgfile "$outputdir/genome.bg"     \
        --topn "$topn"                      \
        --topk "$maxk"                      \
        --oc "$outputdir"                   \
        "$meme_file"                        \
        "$genomefile_temp"                  \
        "$outputdir/promoter_lengths.txt" &
    [ $((n % threads)) -eq 0 ] && wait
done
wait

# Post-process FIMO output - restore ':' in sequence names
for fimofile in "$outputdir/fimohits"/*.txt; do
    sed -i '' 's/__COLON__/:/g' "$fimofile"
done

sed -i '' 's/__COLON__/:/g' "$outputdir/promoter_lengths.txt"
sed -i '' 's/__COLON__/:/g' "$universefile"



rm -f "$genomefile_temp"
print_orange "   └─ Restored sequence names and cleaned up temporary file"

print_elapsed_time "$step_start"

# ==============================================================================
# SUMMARY
# ==============================================================================

print_green "\n╔══════════════════════════════════════════════════════════════════╗"
print_green "║                          SUMMARY                                 ║"
print_green "╚══════════════════════════════════════════════════════════════════╝"
print_white "  Output directory       : "; print_orange "$outputdir"
print_green "────────────────────────────────────────────────────────────────────"
print_white "  Generated files:\n"

for file in promoter_lengths.txt IC.txt binomial_thresholds.txt universe.txt; do
    if [ -f "$outputdir/$file" ]; then
        print_green "    ✓ $file"
    else
        print_red "    ✗ $file (missing)"
    fi
done

print_green "────────────────────────────────────────────────────────────────────"
print_elapsed_time_with_status "$start_time" "success"
print_green "\n✓ PMET indexing completed successfully!\n"
