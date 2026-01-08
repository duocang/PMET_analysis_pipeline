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
print_middle(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    # 获取终端的宽度
    COLUMNS=$(tput cols)
    # 遍历每一行
    while IFS= read -r line; do
        # 计算需要的空格数来居中文本
        padding=$(( (COLUMNS - ${#line}) / 2 ))
        printf "%${padding}s" ''
        printf "${FLUORESCENT_YELLOW}${line}${NC}\n"
    done <<< "$1"
}

script_dir=$(cd -- "$(dirname "$0")" && pwd)
cd "$script_dir"
data_dir="$script_dir/data"
fetch_script="$script_dir/scripts/fetch_tair10.sh"


echo -e "\n\n"
print_middle "The purpose of this script is to                                      \n"
print_middle "  Search motif pairs on the longest isoform of any genomic elemet       "
print_middle "                                                                      \n"



# Give execute permission to all users for the file.
find . -type f \( -name "*.sh" -o -name "*.pl" \) -exec chmod a+x {} \;

########################## 1. Downloading data #######################################
if [[ ! -s "$data_dir/TAIR10.fasta" || ! -s "$data_dir/TAIR10.gff3" ]]; then
    print_fluorescent_yellow "Downloading genome and annotation...\n"
    bash "$fetch_script"
else
    print_green "Genome and annotation are ready!"
fi


start_time=$SECONDS
################################ 2. input parameters ###################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_genomic_element.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=8
res_dir=results/03_genomic_element


# homotypic
overlap="NoOverlap"
utr="Yes"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=no

# data
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme


# genomic element
echo -e "Select the genomic element:\n    1. 3' UTR\n    2. 5' UTR\n    3. mRNA\n    4. CDS\n    5. Exon"
read -p "Enter your choice (1/2/3/4): " choice
case $choice in
    1) genomic_element="three_prime_UTR"; gff3id='Parent=transcript:' ;;
    2) genomic_element="five_prime_UTR"; gff3id='Parent=transcript:' ;;
    3) genomic_element="mRNA"; gff3id='ID=transcript:';;
    4) genomic_element="CDS"; gff3id='Parent=transcript:' ;;
    5) genomic_element="exon"; gff3id='Parent=transcript:' ;;
    *) echo "Invalid choice. Please enter 1, 2, or 3."; exit 1 ;;
esac
print_fluorescent_yellow "Chosen Genomic Element: $genomic_element"
print_fluorescent_yellow "GFF3 ID format: $gff3id"


# output
homotypic_output=$res_dir/01_homotypic

# heterotypic
task=salt_top300
gene_input_file=data/genes/$task.txt
heterotypic_output=$res_dir/02_heterotypic
icthresh=4

# plot output
plot_output=$res_dir/03_plot

mkdir -p $homotypic_output
mkdir -p $heterotypic_output
mkdir -p $plot_output



############################## 3. Running homotypic #################################
print_green "Running homotypic searching...\n"

$HOMOTYPIC               \
    -r $toolDir          \
    -o $homotypic_output \
    -e $genomic_element  \
    -i $gff3id           \
    -k $maxk             \
    -n $topn             \
    -p $length           \
    -v $overlap          \
    -u $utr              \
    -f $fimothresh       \
    -t $threads          \
    -d $delete_temp      \
    $genome              \
    $anno                \
    $meme

# for task in "salt_top300" "random_genes_300" "genes_cell_type_treatment" "gene_cortex_epidermis_pericycle" "heat_top300"; do

#     heterotypic_output=${heterotypic_output}_${task}
#     plot_output=${plot_output}_${task}
#     gene_input_file=data/genes/$task.txt
#     mkdir -p $heterotypic_output
#     mkdir -p $plot_output

#     ############################ 4. Running heterotypic ###############################
#     print_green "\n\nSearching for heterotypic motif hits...\n"

#     # remove genes not present in pre-computed pmet index
#     grep -Ff $homotypic_output/universe.txt $gene_input_file > $heterotypic_output/new_genes_temp.txt


#     $HETEROTYPIC                                     \
#         -d .                                         \
#         -g $heterotypic_output/new_genes_temp.txt    \
#         -i $icthresh                                 \
#         -p $homotypic_output/promoter_lengths.txt    \
#         -b $homotypic_output/binomial_thresholds.txt \
#         -c $homotypic_output/IC.txt                  \
#         -f $homotypic_output/fimohits                \
#         -o $heterotypic_output                       \
#         -t $threads > $heterotypic_output/pmet.log

#     rm $heterotypic_output/new_genes_temp.txt
#     # merge pmet result
#     cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
#     rm $heterotypic_output/temp*.txt

#     end_time=$SECONDS
#     elapsed_time=$((end_time - start_time))
#     days=$((elapsed_time/86400))
#     hours=$(( (elapsed_time%86400)/3600 ))
#     minutes=$(( (elapsed_time%3600)/60 ))
#     seconds=$((elapsed_time%60))
#     print_orange "      Time taken: $days day $hours hour $minutes minute $seconds second\n"


#     print_green "DONE: heterotypic search"

#     ##################################### Heatmap ##################################
#     print_green "\n\nCreating heatmap...\n"

#     Rscript 05_heatmap.R                        \
#         All                                  \
#         $plot_output/heatmap.png             \
#         $heterotypic_output/motif_output.txt \
#         15                                    \
#         3                                    \
#         6                                    \
#         FALSE

#     Rscript 05_heatmap.R                           \
#         Overlap                                 \
#         $plot_output/heatmap_overlap_unique.png \
#         $heterotypic_output/motif_output.txt    \
#         15                                       \
#         3                                       \
#         6                                       \
#         TRUE

#     Rscript 05_heatmap.R                        \
#         Overlap                              \
#         $plot_output/heatmap_overlap.png     \
#         $heterotypic_output/motif_output.txt \
#         15                                    \
#         3                                    \
#         6                                    \
#         FALSE
# done
# # method
# # filename
# # pmet.out
# # topn
# # histgram_ncol
# # histgram_width
# # unique_cmbination