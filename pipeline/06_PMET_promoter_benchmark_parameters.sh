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
################################# Parameters #######################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_promoters_benchmark_parameters.sh
HETEROTYPIC=$toolDir/pmetParallel
chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=8
res_dir=results/06_PMET_promoter_benchmark_parameters

gff3id="gene_id="
# overlap="NoOverlap"
overlap="Yes"
utr="No"
topnRangeStr="5000"
maxkRangeStr="1-2-3-4-5-6-7-8-9"
promlengthRangeStr="50-200-500-1000-1500-2500-5000"
# data
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/ArabidopsisPBM_20140210.meme
# output
homotypic_output=$res_dir/01_homotypic

# Heterotypic
task=genes_cell_type_treatment
gene_file=data/genes/$task.txt
fimothresh=005
icthresh=4

heterotypic_output=$res_dir/02_heterotypic
logDir=results/logs
plot_output=$res_dir/heatmap

mkdir -p $homotypic_output
mkdir -p $heterotypic_output
mkdir -p $logDir
mkdir -p $plot_output

start_time=$SECONDS
################################## Homotypic #######################################
$HOMOTYPIC                 \
    -r $toolDir            \
    -o $homotypic_output   \
    -i $gff3id             \
    -v $overlap            \
    -u $utr                \
    -t $threads            \
    -n $topnRangeStr       \
    -k $maxkRangeStr       \
    -p $promlengthRangeStr \
    $genome                \
    $anno                  \
    $meme
################################## Heterotypic #######################################
# 定义函数以将逗号分隔的字符串转换为数组
string_to_array() {
    local IFS='-' # 设置字段分隔符为逗号
    read -ra arr <<< "$1" # 读取字符串并分割到数组
    echo "${arr[@]}" # 返回数组
}
topnRange=($(string_to_array "$topnRangeStr"))
maxkRange=($(string_to_array "$maxkRangeStr"))
promlengthRange=($(string_to_array "$promlengthRangeStr"))

for task in "genes_cell_type_treatment" "gene_cortex_epidermis_pericycle" "salt_top300" "heat_top300"; do

    print_green "Genes: ${task}.txt"
    gene_file=data/genes/$task.txt

    for promlength in ${promlengthRange[@]}; do
        for maxk in ${maxkRange[@]}; do
            for topn in ${topnRange[@]}; do

                print_orange             "    parameters:"
                print_fluorescent_yellow "        genes     : $gene_file"
                print_fluorescent_yellow "        promlength: $promlength"
                print_fluorescent_yellow "        maxk      : $maxk"
                print_fluorescent_yellow "        topn      : $topn"
                print_fluorescent_yellow "        fimothresh: $fimothresh"
                print_green              "    Searching for heterotypic motifs pairs with ${threads} CPUS..."

                indexOutputTemp=$homotypic_output/fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}
                outputTemp=$heterotypic_output/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}

                mkdir -p $outputTemp

                # remove genes not present in pre-computed pmet index
                grep -Ff $indexOutputTemp/universe.txt $gene_file > $outputTemp/available_genes.txt

                $HETEROTYPIC \
                    -d .                                         \
                    -g $outputTemp/available_genes.txt           \
                    -i $icthresh                                 \
                    -p $indexOutputTemp/promoter_lengths.txt     \
                    -b $indexOutputTemp/binomial_thresholds.txt  \
                    -c $indexOutputTemp/IC.txt                   \
                    -f $indexOutputTemp/fimohits                 \
                    -o $outputTemp                               \
                    -t $threads > $logDir/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}.log 2>&1

                rm  $outputTemp/available_genes.txt
                cat $outputTemp/*.txt > $outputTemp/motif_output.txt
                rm  $outputTemp/temp*.txt
            done
        done
    done
done
end_time=$SECONDS
elapsed_time=$((end_time - start_time))
days=$((elapsed_time/86400))
hours=$(( (elapsed_time%86400)/3600 ))
minutes=$(( (elapsed_time%3600)/60 ))
seconds=$((elapsed_time%60))
print_orange "      Time taken: $days day $hours hour $minutes minute $seconds second\n"
