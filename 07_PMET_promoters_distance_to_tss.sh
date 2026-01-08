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

################################ 1. Downloading data #######################################
if [[ ! -s "$data_dir/TAIR10.fasta" || ! -s "$data_dir/TAIR10.gff3" ]]; then
    print_green "Downloading genome and annotation...\n"
    bash "$fetch_script"
else
    print_green "Genome and annotation are ready!"
fi

################################ 2. input parameters ###################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_promoters_with_distance_tss_fimo_intergrated.sh
HETEROTYPIC=$toolDir/pmetParallel
chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=4
res_dir=results/07_distance_to_tss

gff3id="gene_id="
# overlap="NoOverlap"
overlap="Yes"
utr="Yes"
topn=5000
maxk=5
length=300
fimothresh=0.05
gap=100
promlengthlimit=300
gff3id="gene_id="
delete_temp=no
# data
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme
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
#################################### 3. Running PMET ########################################
for promlengthlimit in 20; do
    for length in 50 ; do
        for gap in 200; do
            ################################ Running homotypic ###################################
            homo_temp=$homotypic_output/plen${length}_gap${gap}_plenmin${promlengthlimit}
            mkdir -p $homo_temp

            print_green "Running fimo...\n"
            $HOMOTYPIC                \
                -r $toolDir           \
                -o $homo_temp         \
                -i $gff3id            \
                -k $maxk              \
                -n $topn              \
                -p $length            \
                -g $gap               \
                -l $promlengthlimit   \
                -v $overlap           \
                -u $utr               \
                -f $fimothresh        \
                -t $threads           \
                -d $delete_temp       \
                $genome               \
                $anno                 \
                $meme
            for task in "genes_cell_type_treatment" "gene_cortex_epidermis_pericycle"; do
                gene_input_file=data/genes/$task.txt

                heterotypic_output=$res_dir/02_heterotypic_${task}
                hetero_temp=$heterotypic_output/plen${length}_gap${gap}_plenmin${promlengthlimit}
                plot_output=$res_dir/plot_${task}/heatmap

                mkdir -p $hetero_temp
                mkdir -p $plot_output

                # ########################## Running heterotypic ##################################
                # print_green "\n\nSearching for heterotypic motif hits..."
                # # remove genes not present in pre-computed pmet index
                # grep -Ff $homo_temp/universe.txt $gene_input_file > $gene_input_file"temp"
                # $HETEROTYPIC                              \
                #     -d .                                  \
                #     -g $gene_input_file"temp"             \
                #     -i $icthresh                          \
                #     -p $homo_temp/promoter_lengths.txt    \
                #     -b $homo_temp/binomial_thresholds.txt \
                #     -c $homo_temp/IC.txt                  \
                #     -f $homo_temp/fimohits                \
                #     -o $hetero_temp                       \
                #     -t $threads > $hetero_temp/pmet.log

                # cat $hetero_temp/*.txt > $hetero_temp/motif_output.txt
                # rm $hetero_temp/temp*.txt
                # rm $gene_input_file"temp"
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