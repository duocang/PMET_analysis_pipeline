#!/bin/bash
set -euo pipefail

# ==================== Setup ====================
# Project root is parent of pipeline directory
script_dir=$(cd -- "$(dirname "$0")/.." && pwd)
cd "$script_dir"
source scripts/lib/print_colors.sh

# ==================== Parameters ====================
# Paths
toolDir=scripts
pmetDir=build
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme

# Homotypic
gff3id="gene_id="
overlap="NoOverlap"
utr="Yes"
topn=5000
maxk=5
length=1000
fimothresh=0.05
delete_temp=no
isPoisson=false

# Heterotypic
task=genes_cell_type_treatment
gene_input_file=data/genes/$task.txt
icthresh=4

# Output
res_dir=results/01_PMET_promoter
homotypic_output=$res_dir/01_homotypic
heterotypic_output=$res_dir/02_heterotypic
plot_output=$res_dir/plot

# Runtime
threads=4

# ==================== Prepare ====================
rm -rf "$homotypic_output" "$heterotypic_output" "$plot_output"
mkdir -p "$homotypic_output" "$heterotypic_output" "$plot_output"
chmod a+x scripts/gff3sort/gff3sort.pl

# Download genome/annotation if missing
if [[ ! -s data/TAIR10.fasta || ! -s data/TAIR10.gff3 ]]; then
    print_green "Downloading genome and annotation..."
    bash scripts/fetch_tair10.sh
else
    print_green "Genome and annotation are ready!"
fi

# ==================== 1. Homotypic ====================
print_green "\n[1/3] Searching for homotypic motif hits..."

bash "$toolDir/PMETindex_promoters_fimo_integrated.sh" \
    -r "$toolDir"          \
    -o "$homotypic_output" \
    -i "$gff3id"           \
    -k "$maxk"             \
    -n "$topn"             \
    -p "$length"           \
    -v "$overlap"          \
    -u "$utr"              \
    -f "$fimothresh"       \
    -t "$threads"          \
    -d "$delete_temp"      \
    -x "$isPoisson"        \
    "$genome" "$anno" "$meme"

rm -rf "$homotypic_output/fimo"

# ==================== 2. Heterotypic ====================
print_green "\n[2/3] Searching for heterotypic motif hits..."

# Filter genes present in index
gene_tmp="${gene_input_file}.tmp"
grep -Ff "$homotypic_output/universe.txt" "$gene_input_file" > "$gene_tmp"

"$pmetDir/pair_parallel" \
    -d .                                           \
    -g "$gene_tmp"                                 \
    -i "$icthresh"                                 \
    -p "$homotypic_output/promoter_lengths.txt"    \
    -b "$homotypic_output/binomial_thresholds.txt" \
    -c "$homotypic_output/IC.txt"                  \
    -f "$homotypic_output/fimohits"                \
    -o "$heterotypic_output"                       \
    -t "$threads" > "$heterotypic_output/pmet.log"

cat "$heterotypic_output"/*.txt > "$heterotypic_output/motif_output.txt"
rm -f "$heterotypic_output"/temp*.txt "$gene_tmp"

# ==================== 3. Heatmap ====================
print_green "\n[3/3] Creating heatmaps..."

Rscript scripts/r/draw_heatmap.R All     "$plot_output/heatmap.png"                "$heterotypic_output/motif_output.txt" 5 3 6 FALSE
Rscript scripts/r/draw_heatmap.R Overlap "$plot_output/heatmap_overlap_unique.png" "$heterotypic_output/motif_output.txt" 5 3 6 TRUE
Rscript scripts/r/draw_heatmap.R Overlap "$plot_output/heatmap_overlap.png"        "$heterotypic_output/motif_output.txt" 5 3 6 FALSE

print_green "\nDone!"