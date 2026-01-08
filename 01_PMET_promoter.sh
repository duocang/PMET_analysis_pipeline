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

script_dir=$(cd -- "$(dirname "$0")" && pwd)
cd "$script_dir"
data_dir="$script_dir/data"
fetch_script="$script_dir/scripts/fetch_tair10.sh"
# Give execute permission to all users for the file.
chmod a+x scripts/cpp_debug_needed/homotypic_promoters.sh
chmod a+x scripts/gff3sort/gff3sort.pl

########################## Downloading data #######################################
if [[ ! -s "$data_dir/TAIR10.fasta" || ! -s "$data_dir/TAIR10.gff3" ]]; then
    print_green "Downloading genome and annotation...\n"
    bash "$fetch_script"
else
    print_green "Genome and annotation are ready!"
fi
########################## Parameters #######################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_promoters_fimo_integrated.sh
HETEROTYPIC=$toolDir/pmetParallel
chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=8
res_dir=results/01_PMET_promoter

# homotypic
gff3id="gene_id="
overlap="NoOverlap"
# overlap="Yes"
utr="Yes"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=no
isPoisson=false
# data
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme
# output
homotypic_output=$res_dir/01_homotypic
# heterotypic
task=genes_cell_type_treatment
gene_input_file=data/genes/$task.txt
heterotypic_output=$res_dir/02_heterotypic
icthresh=4

# plot output
plot_output=$res_dir/plot
mkdir -p $homotypic_output
mkdir -p $heterotypic_output
mkdir -p $plot_output

################################ Running homotypic ###################################
print_green "Searching for homotypic motif hits..."
$HOMOTYPIC               \
    -r $toolDir          \
    -o $homotypic_output \
    -i $gff3id           \
    -k $maxk             \
    -n $topn             \
    -p $length           \
    -v $overlap          \
    -u $utr              \
    -f $fimothresh       \
    -t $threads          \
    -d $delete_temp      \
    -x $isPoisson        \
    $genome              \
    $anno                \
    $meme
rm -rf $homotypic_output/fimo

##################################### Heterotypic ##################################
print_green "\n\nSearching for heterotypic motif hits..."
# remove genes not present in pre-computed pmet index
grep -Ff $homotypic_output/universe.txt $gene_input_file > $gene_input_file"temp"

$HETEROTYPIC                                     \
    -d .                                         \
    -g $gene_input_file"temp"                    \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output                       \
    -t $threads > $heterotypic_output/pmet.log
cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
rm $heterotypic_output/temp*.txt
rm $gene_input_file"temp"


print_green "\n\nSearching for heterotypic motif hits..."
# remove genes not present in pre-computed pmet index
grep -Ff $homotypic_output/universe.txt $gene_input_file > $gene_input_file"temp"

$toolDir/pmet                                    \
    -d .                                         \
    -g $gene_input_file"temp"                    \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output
rm $gene_input_file"temp"

##################################### Heatmap ##################################
print_green "\n\nCreating heatmap...\n"

Rscript 05_heatmap.R                     \
    All                                  \
    $plot_output/heatmap.png             \
    $heterotypic_output/motif_output.txt \
    5                                    \
    3                                    \
    6                                    \
    FALSE

Rscript 05_heatmap.R                        \
    Overlap                                 \
    $plot_output/heatmap_overlap_unique.png \
    $heterotypic_output/motif_output.txt    \
    5                                       \
    3                                       \
    6                                       \
    TRUE

Rscript 05_heatmap.R                     \
    Overlap                              \
    $plot_output/heatmap_overlap.png     \
    $heterotypic_output/motif_output.txt \
    5                                    \
    3                                    \
    6                                    \
    FALSE

# method
# filename
# pmet.out
# topn
# histgram_ncol
# histgram_width
# unique_cmbination