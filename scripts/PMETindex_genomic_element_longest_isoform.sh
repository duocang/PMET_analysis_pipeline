#!/bin/bash
set -e

function usage () {
    cat >&2 <<EOF
USAGE: PMETindexgenome [options] <genome> <gff3> <memefile>

Creates PMET index for Paired Motif Enrichment Test using genome files.
Required arguments:
-r <PMETindex_path>	: Full path of python scripts called from this file. Required.
-i <gff3_identifier> : gene identifier in gff3 file e.g. gene_id=

Optional arguments:
-o <output_directory> : Output directory for results
-n <topn>	: How many top promoter hits to take per motif. Default=5000
-k <max_k>	: Maximum motif hits allowed within each promoter.  Default: 5
-p <promoter_length>	: Length of promoter in bp used to detect motif hits default: 1000
-v <include_overlaps> :  Remove promoter overlaps with gene sequences. AllowOverlap or NoOverlap, Default : AllowOverlap
-u <include_UTR> : Include 5' UTR sequence? Yes or No, default : No
-f <fimo_threshold> : Specify a minimum quality for hits matched by fimo. Default: 0.05

EOF
}

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


# set up defaults
topn=5000
maxk=5
promlength=1000
fimothresh=0.05
element=mRNA
mrnaFull=No
delete=yes
overlap="AllowOverlap"
utr="No"
gff3id='transcript:'
pmetroot="scripts"
threads=4
icthreshold=24

# set up empty variables

indexingOutputDir=
genomefile=
gff3file=
memefile=

# deal with arguments
# if none, exit
if [ $# -eq 0 ]
    then
        echo "No arguments supplied"  >&2
        usage
        exit 1
fi

while getopts ":r:i:o:n:k:p:f:v:u:e:m:t:d:" options; do
    case $options in
        r) pmetroot=$OPTARG;;
        i) gff3id=$OPTARG;;
        o) indexingOutputDir=$OPTARG;;
        n) topn=$OPTARG;;
        k) maxk=$OPTARG;;
        p) promlength=$OPTARG;;
        f) fimothresh=$OPTARG;;
        v) overlap=$OPTARG;;
        u) utr=$OPTARG;;
        e) element=$OPTARG;;
        m) mrnaFull=$OPTARG;;
        t) threads=$OPTARG;;
        d) delete=$OPTARG;;
        \?) print_red  "Invalid option: -$OPTARG" >&2
        exit 1;;
        :)  print_red "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done


shift $((OPTIND - 1))
genomefile=$1
gff3file=$2
memefile=$3
universefile=$indexingOutputDir/universe.txt
bedfile=$indexingOutputDir/${element}.bed


print_white "Genome file                  : "; print_orange $genomefile
print_white "Annotation file              : "; print_orange $gff3file
print_white "Motif meme file              : "; print_orange $memefile

print_white "PMET index path              : "; print_orange "$pmetroot"
print_white "GFF3 identifier              : "; print_orange "$gff3id"
print_white "Genomic element              : "; print_orange "$element"
print_white "Output directory             : "; print_orange "$indexingOutputDir"
print_white "Top n promoters              : "; print_orange "$topn"  # Default to 5000 if not set
print_white "Top k motif hits             : "; print_orange "$maxk"     # Default to 5 if not set
print_white "Length of promoter           : "; print_orange "$promlength"  # Default to 1000 if not set
print_white "Fimo threshold               : "; print_orange "$fimothresh"
print_white "Promoter overlap handling    : "; print_orange "$overlap"
print_white "Include 5' UTR               : "; print_orange "$utr"
print_white "Number of threads            : "; print_orange "$threads"

mkdir -p $indexingOutputDir

start=$SECONDS

print_green "Preparing data for FIMO and PMET index..."

# -------------------------------------------------------------------------------------------
# 1. sort annotaion by gene coordinates
print_fluorescent_yellow "     1. Sorting annotation by gene coordinates"
$pmetroot/gff3sort/gff3sort.pl $gff3file > $indexingOutputDir/sorted.gff3

# -------------------------------------------------------------------------------------------
# 2. extract gene line from annoitation
print_fluorescent_yellow "     2. Extracting gene line from annoitation"
# grep -P '\mRNA\t' $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
if [[ "$(uname)" == "Linux" ]]; then
    grep -P "\t${element}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
elif [[ "$(uname)" == "Darwin" ]]; then
    grep    "\t${element}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
else
    print_red "Unsupported operating system."
fi

# -------------------------------------------------------------------------------------------
# 3. extract chromosome , start, end, gene ('gene_id' for input) ...
print_fluorescent_yellow "     3. Extracting chromosome, start, end, gene"
python3 $pmetroot/parse_mRNAlines.py $gff3id $indexingOutputDir/genelines.gff3 $bedfile

# -------------------------------------------------------------------------------------------
# 4. filter invalid genes: start should be smaller than end
invalidRows=$(awk '$2 >= $3' $bedfile)
if [[ -n "$invalidRows" ]]; then
    echo "$invalidRows" > $indexingOutputDir/invalid_lines.bed
fi
print_fluorescent_yellow "     4. Extracting genes coordinates: start should be smaller than end (genelines.bed)"
# 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
# In the BED file format, the start position is always less than the end position for both positive (+) and negative (-) chains.
# 起始和终止位置是指定基因上的物理位置，而不是表达或翻译的方向。
# start and end positions specify the physical location of the gene, rather than the direction of expression or translation.
# starting site < stopped site in bed file
awk '$2 <  $3' $bedfile > temp.bed && mv temp.bed $bedfile

# -------------------------------------------------------------------------------------------
# 5. convert isoform name to it gene
print_fluorescent_yellow "     5. Converting isoform name to its gene (cleaned.bed)"
awk -v OFS="\t" '{split($4, arr, "."); $4 = arr[1]; print $0}' $bedfile > $indexingOutputDir/cleaned.bed

# -------------------------------------------------------------------------------------------
# 6. find longest genes (it was isoforms before removing .x)
print_fluorescent_yellow "     6. Finding longest isoform"
awk '
{
    len = $3 - $2;
    if (!($4 in maxLen) || len > maxLen[$4]) {
        maxLen[$4] = len;  # 更新最大长度
        maxLine[$4] = $0;  # 保留整行内容
    }
}
END {
    for (id in maxLine) {
        printf "%s\n", maxLine[id];
    }
}' $indexingOutputDir/cleaned.bed > $indexingOutputDir/cleaned_longest.bed

sort -k4,4 $indexingOutputDir/cleaned_longest.bed > $bedfile
rm -rf $indexingOutputDir/cleaned_longest.bed
rm -rf $indexingOutputDir/cleaned.bed


if [ "$element" = "mRNA" ] && [ "$mrnaFull" = "No" ]; then
    print_fluorescent_yellow "        (only for mRNA) Removing overlap with 3' UTR and 5' UTR"
    for i in 3 5; do
        if [ "$i" -eq 3 ]; then
            item="three_prime_UTR"
            gff3id="Parent=transcript:"
        else
            item="five_prime_UTR"
            gff3id="Parent=transcript:"
        fi
        # 2. extract gene line from annoitation
        if [[ "$(uname)" == "Linux" ]]; then
            grep -P "\t${item}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/${item}.gff3
        elif [[ "$(uname)" == "Darwin" ]]; then
            grep    "\t${item}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/${item}.gff3
        fi
        # 3. extract chromosome , start, end, gene ('gene_id' for input) ...
        python3 $pmetroot/parse_mRNAlines.py $gff3id $indexingOutputDir/${item}.gff3 $indexingOutputDir/${item}.bed
        # 4. filter invalid genes: start should be smaller than end
        awk '$2 <  $3' $indexingOutputDir/${item}.bed > temp.bed && mv temp.bed $indexingOutputDir/${item}.bed
        # 5. convert isoform name to it gene
        awk -v OFS="\t" '{split($4, arr, "."); $4 = arr[1]; print $0}' $indexingOutputDir/${item}.bed > $indexingOutputDir/${item}_cleaned.bed
        # 6. find longest genes (it was isoforms before removing .x)
        awk '
        {
            len = $3 - $2;
            if (!($4 in maxLen) || len > maxLen[$4]) {
                maxLen[$4] = len;  # 更新最大长度
                maxLine[$4] = $0;  # 保留整行内容
            }
        }
        END {
            for (id in maxLine) {
                printf "%s\n", maxLine[id];
            }
        }' $indexingOutputDir/${item}_cleaned.bed > $indexingOutputDir/${item}_cleaned_longest.bed
        sort -k4,4 $indexingOutputDir/${item}_cleaned_longest.bed > $indexingOutputDir/${item}.bed
        rm -rf $indexingOutputDir/${item}_cleaned_longest.bed
        rm -rf $indexingOutputDir/${item}_cleaned.bed
        rm -rf $indexingOutputDir/${item}.gff3

        bedtools subtract -a  $bedfile -b $indexingOutputDir/${item}.bed > $indexingOutputDir/non_overlapping.bed
    done

    cp $bedfile $indexingOutputDir/with_overlapping.bed

    bedtools subtract                             \
        -a $bedfile                               \
        -b $indexingOutputDir/three_prime_UTR.bed \
        > $indexingOutputDir/temp.bed
    bedtools subtract                             \
        -a $indexingOutputDir/temp.bed            \
        -b $indexingOutputDir/five_prime_UTR.bed  \
        > $indexingOutputDir/non_overlapping.bed

    rm -rf $indexingOutputDir/temp.bed
    rm -rf $bedfile
    mv $indexingOutputDir/non_overlapping.bed $bedfile

    print_fluorescent_yellow "        (only for mRNA) Makiing fragments of mRNA to be an independent mRNA"
    # mRNA can be divided by removing 3' and 5'
    cut -f4 $bedfile | sort | uniq -d > $indexingOutputDir/id_duplicated.txt
    # add '__' and '__[number]' to make mRNA fragments unique: AT5g139303 and __AT5g139303__2
    awk '{
        if (seen[$4]++ == 0) {
            # 如果是第一次看到这个基因名称，不做改变
            print $0;
        } else {
            # 对于重复的基因名称，在名称前加 "__"，后面加 "__" 和一个序号
            print $1 "\t" $2 "\t" $3 "\t" "__" $4 "__" seen[$4] "\t" $5 "\t" $6;
        }
    }' "$bedfile" >  $indexingOutputDir/modified_bedfile.bed
    # filter mRNA with too short sequences
    awk '$3 - $2 >= 30' $indexingOutputDir/modified_bedfile.bed > $indexingOutputDir/filtered_bedfile.bed

    sort -k4,4 $indexingOutputDir/filtered_bedfile.bed \
        > $indexingOutputDir/sorted_filtered_bedfile.bed

    rm -rf $bedfile
    mv $indexingOutputDir/sorted_filtered_bedfile.bed $bedfile
fi



# -------------------------------------------------------------------------------------------
# 7. list of all genes found
print_fluorescent_yellow "     7. Extracting all genes found (universe.txt)"
cut -f 4 $bedfile > $universefile

# -------------------------------------------------------------------------------------------
# 8. promoter lenfths from promoters.bed
print_fluorescent_yellow "     8. Promoter lengths from genelines.bed"
awk '{print $4 "\t" ($3 - $2)}' $bedfile      \
    > $indexingOutputDir/promoter_lengths.txt

# -------------------------------------------------------------------------------------------
# 9. strip the potential FASTA line breaks. creates genome_stripped.fa
print_fluorescent_yellow "     9. Removing potential FASTA line breaks (genome_stripped.fa)"
awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
    { printf "%s",$0;} \
    END { print ""; }'  $genomefile > $indexingOutputDir/genome_stripped.fa


# -------------------------------------------------------------------------------------------
# 10. create promoters fasta
print_fluorescent_yellow "    10. Creating mRNA FASTA file (promoter_rought.fa)";
bedtools getfasta \
        -fi  $indexingOutputDir/genome_stripped.fa \
        -bed $bedfile                              \
        -fo  $indexingOutputDir/promoter_rought.fa \
        -name

# -------------------------------------------------------------------------------------------
# 11. replace the id of each seq with gene names
print_fluorescent_yellow "    11. Replacing the id of each sequences' with gene names (promoter.fa)"
sed 's/::.*//g' $indexingOutputDir/promoter_rought.fa > $indexingOutputDir/promoter.fa

# check if any duplicated id
grep "^>" $indexingOutputDir/promoter.fa > $indexingOutputDir/ids.txt
sort $indexingOutputDir/ids.txt | uniq -d > $indexingOutputDir/duplicate_ids.txt
if [ ! -s $indexingOutputDir/duplicate_ids.txt ]; then
    rm -rf $indexingOutputDir/ids.txt
    rm -rf $indexingOutputDir/duplicate_ids.txt
fi

# -------------------------------------------------------------------------------------------
# 12. promoter.bg from promoter.fa
print_fluorescent_yellow "    12. fasta-get-markov estimates a Markov model from promoter.fa. (promoter.bg)"
fasta-get-markov $indexingOutputDir/promoter.fa > $indexingOutputDir/promoter.bg

# -------------------------------------------------------------------------------------------
# 13. IC.txt
print_fluorescent_yellow "    13. Generating information content (IC.txt)"
[ ! -d $indexingOutputDir/memefiles ] && mkdir $indexingOutputDir/memefiles
python3 $pmetroot/parse_memefile.py $memefile $indexingOutputDir/memefiles/
python3 $pmetroot/calculateICfrommeme_IC_to_csv.py \
    $indexingOutputDir/memefiles/                  \
    $indexingOutputDir/IC.txt
rm -rf $indexingOutputDir/memefiles/*

# -------------------------------------------------------------------------------------------
# 14. individual motif files from user's meme file
print_fluorescent_yellow "    14. Spliting motifs into individual meme files (folder memefiles)"
python3 $pmetroot/parse_memefile_batches.py $memefile $indexingOutputDir/memefiles/ $threads

# -------------------------------- Run fimo and pmetindex --------------------------
[ ! -d $indexingOutputDir/fimohits ] && mkdir $indexingOutputDir/fimohits

print_green "Running FIMO and PMET index..."
runFimoIndexing () {
    memefile=$1
    indexingOutputDir=$2
    fimothresh=$3
    pmetroot=$4
    maxk=$5
    topn=$6
    filename=`basename $memefile .txt`

    $pmetroot/fimo                              \
        --no-qvalue                             \
        --text                                  \
        --thresh $fimothresh                    \
        --verbosity 1                           \
        --bgfile $indexingOutputDir/promoter.bg \
        --topn $topn                            \
        --topk $maxk                            \
        --oc $indexingOutputDir/fimohits        \
        $memefile                               \
        $indexingOutputDir/promoter.fa          \
        $indexingOutputDir/promoter_lengths.txt
}
export -f runFimoIndexing


nummotifs=$(grep -c '^MOTIF' "$memefile")
print_orange "    $nummotifs motifs found"

find $indexingOutputDir/memefiles -name \*.txt \
    | parallel --progress --jobs=$threads \
        "runFimoIndexing {} $indexingOutputDir $fimothresh $pmetroot $maxk $topn"

mv $indexingOutputDir/fimohits/binomial_thresholds.txt $indexingOutputDir/



if [ "$element" = "mRNA" ] && [ "$mrnaFull" = "No" ]; then
    print_fluorescent_yellow "        (only for mRNA) Adding up lengths of mRNA with multiple fragments"
    # remove '__' and '__[numbe]' of a gene (mRAN frangments): __AT5g139303__2
    awk '{
        sub(/^__/, "", $1);
        sub(/__[0-9]+$/, "", $1);
        print $1 "\t" $2;
    }' $indexingOutputDir/promoter_lengths.txt > $indexingOutputDir/promoter_lengths_modified.txt
    # add up lengths and keep only one genes for multiple mRAN fragments
    awk '{
    sum[$1] += $2; # 对于每个基因，累加第二列的值
    }
    END {
        for (gene in sum) {
            print gene "\t" sum[gene]; # 打印每个基因及其对应的总和
        }
    }' $indexingOutputDir/promoter_lengths_modified.txt > $indexingOutputDir/promoter_lengths_summed.txt

    rm -rf $indexingOutputDir/promoter_lengths.txt
    mv $indexingOutputDir/promoter_lengths_summed.txt $indexingOutputDir/promoter_lengths.txt

    rm -rf $indexingOutputDir/promoter_lengths_summed.txt
    rm -rf $indexingOutputDir/promoter_lengths_modified.txt

    print_fluorescent_yellow "        (only for mRNA) Merging and filtering fimo hits out of multiple mRNA fragments"
    Rscript $pmetroot/parse_mRNA_multiple_fragments.r \
        $indexingOutputDir/binomial_thresholds.txt    \
        $indexingOutputDir/fimohits                   \
        $indexingOutputDir/fimohits_

    rm -rf $indexingOutputDir/fimohits
    mv $indexingOutputDir/fimohits_ $indexingOutputDir/fimohits
fi


# -------------------------------------------------------------------------------------------
# Deleting unnecessary files
if [[ $delete == "yes" || $delete == "YES" || $delete == "Y" || $delete == "y" ]]; then
    print_green "Deleting unnecessary files...\n\n"
    # rm -rf $indexingOutputDir/IC.txt
    # rm -rf $indexingOutputDir/binomial_thresholds.txt
    rm -rf $indexingOutputDir/filtered_bedfile.bed
    # rm -rf $indexingOutputDir/fimohits
    rm -rf $indexingOutputDir/fimohits_
    rm -rf $indexingOutputDir/five_prime_UTR.bed
    rm -rf $indexingOutputDir/genelines.gff3
    rm -rf $indexingOutputDir/genome_stripped.fa
    rm -rf $indexingOutputDir/genome_stripped.fa.fai
    rm -rf $indexingOutputDir/id_duplicated.txt
    rm -rf $indexingOutputDir/mRNA.bed
    # rm -rf $indexingOutputDir/memefiles
    rm -rf $indexingOutputDir/modified_bedfile.bed
    rm -rf $indexingOutputDir/promoter.bg
    rm -rf $indexingOutputDir/promoter.fa
    # rm -rf $indexingOutputDir/promoter_lengths.txt
    rm -rf $indexingOutputDir/promoter_rought.fa
    rm -rf $indexingOutputDir/sorted.gff3
    rm -rf $indexingOutputDir/three_prime_UTR.bed
    # rm -rf $indexingOutputDir/universe.txt
    rm -rf $indexingOutputDir/with_overlapping.bed
fi

# -------------------------------------------------------------------------------------------
# Checking results
print_green "Checking results...\n\n"
# 计算 $indexingOutputDir/fimohits 目录下 .txt 文件的数量
# Count the number of .txt files in the $indexingOutputDir/fimohits directory
file_count=$(find "$indexingOutputDir/fimohits" -maxdepth 1 -type f -name "*.txt" | wc -l)

# 检查文件数量是否等于 meotif的数量 （$nummotifs）
# Check if the number of files equals the number of meotifs ($nummotifs)
if [ "$file_count" -eq "$nummotifs" ]; then
    end=$SECONDS
    elapsed_time=$((end - start))
    days=$((elapsed_time/86400))
    hours=$(( (elapsed_time%86400)/3600 ))
    minutes=$(( (elapsed_time%3600)/60 ))
    seconds=$((elapsed_time%60))
    print_orange "      Time take: $days day $hours hour $minutes minute $seconds seconds"

    print_green "DONE: homotypic search"
else
    print_red "\nError: there are $file_count fimohits files, it should be $nummotifs."
fi

# # next stage needs the following inputs

# #   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
# #   bimnomial_thresholds.txt    made by PMETindex
# #   IC.txt                      made by calculateICfrommeme.py from meme file
# #   gene input file             supplied by user
