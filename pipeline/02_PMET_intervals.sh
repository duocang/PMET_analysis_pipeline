#!/bin/bash

# Project root is parent of pipeline directory
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$script_dir"
print_colors_lib="$script_dir/scripts/lib/print_colors.sh"

if [ ! -f "$print_colors_lib" ]; then
    echo "ERROR: Missing print_colors library at $print_colors_lib" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$print_colors_lib"

##################################### 1. User input #######################################
res_dir=results/02_PMET_intervals

genome=data/homotypic_intervals/intervals.fa
meme=data/homotypic_intervals/motif_more.meme
gene_input_file=data/homotypic_intervals/intervals.txt
# genome=data_temp/1.fa
# meme=data_temp/NKX6_test.meme
# gene_input_file=data_temp/1.txt

##################################### 2. Parameters #######################################
threads=1

toolDir=scripts

# homotypic
gff3id="gene_id="
overlap="NoOverlap"
utr="No"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=yes
icthresh=4

################################### 3. mkdir and path #####################################
rm -rf $res_dir
mkdir -p $res_dir
# output
homotypic_output=$res_dir/01_homotypic
heterotypic_output=$res_dir/02_heterotypic
mkdir -p $homotypic_output
mkdir -p $heterotypic_output

################################# 4. chmod execute permission ############################
# Give execute permission to all users for the file.
# tool path
HOMOTYPIC=$toolDir/cpp_debug_needed/homotypic_intervals.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x scripts/cpp_debug_needed/homotypic_intervals.sh
chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

find ../../01_PMETDEV-code/debug/scripts/ -type f \( -name "*.sh" -o -name "*.perl" \) -exec chmod a+x {} \;

##################################### 5. Homotypic #####################################
print_green "Running homotypic...\n"
echo "Homotypic output: $homotypic_output"
$HOMOTYPIC               \
    -r $toolDir          \
    -o $homotypic_output \
    -k $maxk             \
    -n $topn             \
    -f $fimothresh       \
    -t $threads          \
    $genome              \
    $meme

################################### 6. Heterotypic #####################################
# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' homotypic_output/promoter_lengths.txt > homotypic_output/temp_genes_list.txt
# cat homotypic_output/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm homotypic_output/temp_genes_list.txt

print_green "Searching for heterotypic motif hits..."
echo "Heterotypic output: $heterotypic_output"
$HETEROTYPIC                                     \
    -d .                                         \
    -g $gene_input_file                          \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output                       \
    -t $threads

cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
rm $heterotypic_output/temp*.txt

#################################### Heatmap ##################################
Rscript scripts/r/draw_heatmap.R                     \
    Overlap                              \
    $heterotypic_output/heatmap.png      \
    $heterotypic_output/motif_output.txt \
    5                                    \
    3                                    \
    6                                    \
    FALSE
