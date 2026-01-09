#!/bin/bash

function error_exit() {
    echo "ERROR: $1" >&2
    usage
    exit 1
}

print_red(){
    RED='\033[0;31m'
    NC='\033[0m' # No Color
    printf "${RED}$1${NC}\n"
}

print_green(){
    GREEN='\033[0;32m'
    NC='\033[0m' # No Color
    printf "${GREEN}$1${NC}\n"
}

print_orange(){
    ORANGE='\033[0;33m'
    NC='\033[0m' # No Color
    printf "${ORANGE}$1${NC}\n"
}

print_fluorescent_yellow(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    printf "${FLUORESCENT_YELLOW}$1${NC}\n"
}

print_white(){
    WHITE='\033[1;37m'
    NC='\033[0m' # No Color
    printf "${WHITE}$1${NC}"
}

script_dir=$(cd -- "$(dirname "$0")/.." && pwd)
cd "$script_dir"

output=results/08_heterotypic_promoters_single_CPU
indexoutput=data/homotypic_promoters
gene_input_file=data/gene.txt



print_green "Searching for heterotypic motif hits with single CPU..."
outputTemp=$output/single
mkdir -p $outputTemp
# remove genes not present in pre-computed pmet index
grep -Ff $indexoutput/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmet                                \
    -d .                                    \
    -g $gene_input_file"temp"               \
    -i 4                                    \
    -p $indexoutput/promoter_lengths.txt    \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt                  \
    -f $indexoutput/fimohits                \
    -o $output/single

# rm $gene_input_file"temp"

Rscript scripts/r/draw_heatmap.R         \
    Overlap                  \
    $output/heatmap.png      \
    $output/motif_output.txt



##################################### pmetParallel ##################################
print_green "Searching for heterotypic motif hits with single CPU..."
outputTemp=$output/parallel
mkdir -p $outputTemp

# remove genes not present in pre-computed pmet index
grep -Ff $indexoutput/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmetParallel                        \
    -d .                                    \
    -g $gene_input_file"temp"               \
    -i 4                                    \
    -p $indexoutput/promoter_lengths.txt    \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt                  \
    -f $indexoutput/fimohits                \
    -t 2                                    \
    -o $outputTemp

cat $outputTemp/*.txt > $outputTemp/motif_output.txt
rm $outputTemp/temp*.txt
rm $gene_input_file"temp"



Rscript scripts/r/draw_heatmap.R             \
    Overlap                      \
    $outputTemp/heatmap.png      \
    $outputTemp/motif_output.txt
