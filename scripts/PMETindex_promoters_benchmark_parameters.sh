#!/bin/bash
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
print_light_blue() {
    LIGHT_BLUE='\033[1;34m'
    NC='\033[0m' # No Color
    printf "${LIGHT_BLUE}$1${NC}\n"
}

# set up defaults
# topn=5000
# maxk=5
# promlength=1000
fimothresh=0.05
overlap="AllowOverlap"
utr="No"
gff3id='gene_id'
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

while getopts ":r:i:o:n:k:p:f:v:u:t:" options; do
    case $options in
        r) pmetroot=$OPTARG;;
        i) gff3id=$OPTARG;;
        o) outputDir=$OPTARG;;
        n) topn=$OPTARG;;
        k) maxk=$OPTARG;;
        p) promlength=$OPTARG;;
        f) fimothresh=$OPTARG;;
        v) overlap=$OPTARG;;
        u) utr=$OPTARG;;
        t) threads=$OPTARG;;
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

if [ true ]; then
    # 定义函数以将逗号分隔的字符串转换为数组 Define function to convert comma separated string to array
    string_to_array() {
        local IFS='-' # 设置字段分隔符为逗号
        read -ra arr <<< "$1" # 读取字符串并分割到数组
        echo "${arr[@]}" # 返回数组
    }
    # 检查数组中的元素是否都是数字 Check if the elements in the array are all numbers
    check_all_numbers() {
        local -a array=("$@") # 从所有传入的参数中重建数组
        unset array[${#array[@]}-1] # 移除最后一个元素（参数名称）
        local param_name="${!#}" # 获取最后一个参数作为参数名称

        for item in "${array[@]}"; do
            if ! [[ "$item" =~ ^[0-9]+$ ]]; then
                print_red "Error: '$item' in $param_name is not a number." >&2
                exit 1
            fi
        done
    }

    topnRange=($(string_to_array "$topn"))
    maxkRange=($(string_to_array "$maxk"))
    promlengthRange=($(string_to_array "$promlength"))

    # 检查每个数组 Check each array
    check_all_numbers "${topnRange[@]}" "topn"
    check_all_numbers "${maxkRange[@]}" "maxk"
    check_all_numbers "${promlengthRange[@]}" "promlength"


    # 检查fimothresh是否为小于1的小数 Check if fimothresh is a decimal less than 1
    if ! [[ "$fimothresh" =~ ^0?\.[0-9]+$ ]] || (( $(echo "$fimothresh >= 1" | bc -l) )); then
        print_red "Error: fimo threshold must be a decimal number less than 1." >&2
        exit 1
    fi

    # 检查输入是否为大于等于1的整数
    if ! [[ "$threads" =~ ^[0-9]+$ ]] || (( threads < 1 )); then
        print_red "Error: The number of threads must be an integer greater than or equal to 1." >&2
        exit 1
    fi

    overlap_lower=$(echo "$overlap" | tr '[:upper:]' '[:lower:]')
    # 检查输入是否为指定的选项之一
    if [[ "$overlap_lower" != "y" && "$overlap_lower" != "n" && "$overlap_lower" != "yes" && "$overlap_lower" != "no" && "$overlap_lower" != "nooverlap" ]]; then
        print_red "Error: Invalid overlap input. Please enter 'y', 'n', 'yes', 'no', or 'NoOverlap'." >&2
        exit 1
    fi

    utr_lower=$(echo "$utr" | tr '[:upper:]' '[:lower:]')
    # 检查输入是否为指定的选项之一
    if [[ "$utr_lower" != "y" && "$utr_lower" != "n" && "$utr_lower" != "yes" && "$utr_lower" != "no" ]]; then
        print_red "Error: Invalid utr input. Please enter 'y', 'n', 'yes', 'no', or 'NoOverlap'." >&2
        exit 1
    fi

fi


print_white "Genome file                  : "; print_orange "$genomefile"
print_white "Annotation file              : "; print_orange "$gff3file"
print_white "Motif meme file              : "; print_orange "$memefile"
print_white "PMET index path              : "; print_orange "$pmetroot"
print_white "GFF3 identifier              : "; print_orange "$gff3id"
print_white "Output directory             : "; print_orange "$outputDir"
print_white "Top n promoters              : "; print_orange "${topn:-5000}"  # Default to 5000 if not set
print_white "Top k motif hits             : "; print_orange "${maxk:-5}"     # Default to 5 if not set
print_white "Length of promoter           : "; print_orange "${promlength:-1000}"  # Default to 1000 if not set
print_white "Fimo threshold               : "; print_orange "$fimothresh"
print_white "Promoter overlap handling    : "; print_orange "$overlap"
print_white "Include 5' UTR               : "; print_orange "$utr"
print_white "Number of threads            : "; print_orange "$threads"

start=$SECONDS

print_green "\nPreparing data for FIMO and PMET index..."

for promlength in ${promlengthRange[@]}; do

    print_green "Preparing data for FIMO and PMET index, promoter legnth: $promlength..."

    indexingOutputDir=$outputDir/LEN${promlength}_FIMO${fimothresh//./}
    universefile=$indexingOutputDir/universe.txt
    bedfile=$indexingOutputDir/genelines.bed

    mkdir -p $indexingOutputDir/fimohits
    # -------------------------------------------------------------------------------------------
    # 1. sort annotaion by gene coordinates
    print_fluorescent_yellow "     1.  Sorting annotation by gene coordinates"
    chmod a+x $pmetroot/gff3sort/gff3sort.pl
    $pmetroot/gff3sort/gff3sort.pl $gff3file > $indexingOutputDir/sorted.gff3

    # -------------------------------------------------------------------------------------------
    # 2. extract gene line from annoitation
    print_fluorescent_yellow "     2.  Extracting gene line from annoitation"
    # grep -P '\tgene\t' $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
    if [[ "$(uname)" == "Linux" ]]; then
        grep -P '\tgene\t' $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
    elif [[ "$(uname)" == "Darwin" ]]; then
        grep '\tgene\t' $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
    else
        print_red "Unsupported operating system."
    fi

    # -------------------------------------------------------------------------------------------
    # 3. extract chromosome , start, end, gene ('gene_id' for input)
    print_fluorescent_yellow "     3.  Extracting chromosome, start, end, gene"

    # 使用grep查找字符串 check if gene_id is present
    if grep -q "$gff3id" "$indexingOutputDir/genelines.gff3"; then
        python3 $pmetroot/parse_genelines.py $gff3id $indexingOutputDir/genelines.gff3 $bedfile
    else
        gff3id='ID='
        python3 $pmetroot/parse_genelines.py $gff3id $indexingOutputDir/genelines.gff3 $bedfile
    fi

    # -------------------------------------------------------------------------------------------
    # 4. filter invalid genes: start should be smaller than end
    print_fluorescent_yellow "     4.  Filter invalid coordinates: start > end"
    touch $indexingOutputDir/invalid_genelines.txt
    awk '$2 >= $3' $bedfile > $indexingOutputDir/invalid_gff3_lines.txt
    awk '$2 <  $3' $bedfile > temp.bed && mv temp.bed $bedfile
    # 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
    # In the BED file format, the start position is always less than the end position for both positive (+) and negative (-) chains.
    # 起始和终止位置是指定基因上的物理位置，而不是表达或翻译的方向。
    # start and end positions specify the physical location of the gene, rather than the direction of expression or translation.



    # -------------------------------------------------------------------------------------------
    print_fluorescent_yellow "         Calculate lenght of space to TSS (length_to_tss.txt)"
    # get length of region before TSS of a gene
    # 初始化变量
    prev_end=0
    prev_chr=""
    next_start=0
    current_line=""

    # 读取BED文件
    while read -r line; do
        if [[ -n $current_line ]]; then
            # 分割当前行和下一行
            read -r chr start end gene score strand <<< "$current_line"
            read -r next_chr next_start next_end next_gene next_score next_strand <<< "$line"
            # 检查染色体是否变化
            if [[ $chr != $prev_chr ]]; then
                prev_end=0
                prev_chr=$chr
            fi
            # 计算与上一个基因末尾的距离
            if [[ $strand == "+" ]]; then
                distance=$((start - prev_end))
            else
                if [[ $chr == $next_chr ]]; then
                    distance=$((next_start - end))
                else
                    distance=0
                fi
            fi
            echo "$gene  $distance" >> $indexingOutputDir/length_to_tss.txt
            # 更新前一个基因的结束位置
            prev_end=$end
        fi
        current_line=$line
    done < $bedfile

    # 处理文件的最后一行
    if [[ -n $current_line ]]; then
        read -r chr start end gene score strand <<< "$current_line"
        if [[ $strand == "+" ]]; then
            distance=$((start - prev_end))
        else
            distance=0
        fi
        echo "$gene  $distance" >> $indexingOutputDir/length_to_tss.txt
    fi
    # draw histogram
    Rscript $pmetroot/histgram_len_to_tss.R $indexingOutputDir/length_to_tss.txt

    # -------------------------------------------------------------------------------------------
    # 5. list of all genes found
    print_fluorescent_yellow "\n     5.  Extracting genes names: complete list of all genes found (universe.txt)"
    cut -f 4 $bedfile > $universefile

    # -------------------------------------------------------------------------------------------
    # 6. strip the potential FASTA line breaks. creates genome_stripped.fa
    print_fluorescent_yellow "     6.  Removing potential FASTA line breaks (genome_stripped.fa)"
    awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
        { printf "%s",$0;} \
        END { print ""; }'  $genomefile > $indexingOutputDir/genome_stripped.fa
    # python3 $pmetroot/strip_newlines.py $genomefile $indexingOutputDir/genome_stripped_py.fa


    # -------------------------------------------------------------------------------------------
    # 7. create the .genome file which contains coordinates for each chromosome start
    print_fluorescent_yellow "     7.  Listing chromosome start coordinates (bedgenome.genome)"
    samtools faidx $indexingOutputDir/genome_stripped.fa
    cut -f 1-2 $indexingOutputDir/genome_stripped.fa.fai > $indexingOutputDir/bedgenome.genome

    # -------------------------------------------------------------------------------------------
    # 8. create promoters' coordinates from annotation
    print_fluorescent_yellow "     8.  Creating promoters' coordinates from annotation (promoters.bed)"
    # 在bedtools中，flank是一个命令行工具，用于在BED格式的基因组坐标文件中对每个区域进行扩展或缩短。
    # In bedtools, flank is a command-line tool used to extend or shorten each region in a BED format genomic coordinate file.
    # 当遇到负链（negative strand）时，在区域的右侧进行扩展或缩短，而不是左侧。
    # When a negative strand is encountered, it is expanded or shortened on the right side of the region, not the left.
    bedtools flank                             \
        -l $promlength                         \
        -r 0 -s -i $bedfile                    \
        -g $indexingOutputDir/bedgenome.genome \
        > $indexingOutputDir/promoters_not_sorted.bed
    # Sort by starting coordinate
    sortBed -i $indexingOutputDir/promoters_not_sorted.bed > $indexingOutputDir/promoters.bed
    rm -rf $indexingOutputDir/promoters_not_sorted.bed

    # # -------------------------------------------------------------------------------------------
    # print_fluorescent_yellow "     8.1 Remove promoters with less than 20 base pairs"
    # # remove promoter length < 20
    # awk '($3 - $2) >= 10' $indexingOutputDir/promoters.bed > $indexingOutputDir/promoters_.bed
    # mv $indexingOutputDir/promoters_.bed $indexingOutputDir/promoters.bed
    # awk '($3 - $2) <  10' $indexingOutputDir/promoters.bed > $indexingOutputDir/8_promoters_less_20.bed

    # -------------------------------------------------------------------------------------------
    # 9. remove overlapping promoter chunks
    if [[ $overlap == 'NoOverlap' || $overlap == "no" || $overlap == "NO" || $overlap == "No" || $overlap == "N" || $overlap == "n" ]]; then
        print_fluorescent_yellow "     9.  Removing overlapping promoter chunks (promoters.bed)"
        sleep 0.1
        bedtools subtract                       \
            -a $indexingOutputDir/promoters.bed \
            -b $bedfile                         \
            > $indexingOutputDir/promoters2.bed
        mv $indexingOutputDir/promoters2.bed $indexingOutputDir/promoters.bed
    else
        print_fluorescent_yellow "     9.  (skipped) Removing overlapping promoter chunks (promoters.bed)"
    fi

    # # -------------------------------------------------------------------------------------------
    # # remove promoter length < 20
    # print_fluorescent_yellow "     9.1 Remove promoters with less than 20 base pairs"
    # awk '($3 - $2) <  10' $indexingOutputDir/promoters.bed > $indexingOutputDir/9_promoters_less_20.bed
    # awk '($3 - $2) >= 10' $indexingOutputDir/promoters.bed > $indexingOutputDir/promoters_.bed
    # mv $indexingOutputDir/promoters_.bed $indexingOutputDir/promoters.bed

    # -------------------------------------------------------------------------------------------
    # 10. check split promoters. if so, keep the bit closer to the TSS
    print_fluorescent_yellow "    10.  Checking split promoter (if so):  keep the bit closer to the TSS (promoters.bed)"
    python3 $pmetroot/assess_integrity.py $indexingOutputDir/promoters.bed

    # TODO: in some case, annotation file sorted by coordinates can have one gene (several lines) splits by other gene (several lines)
    # -------------------------------------------------------------------------------------------
    # 11. add 5' UTR
    if [[ $utr == "yes" || $utr == "YES" || $utr == "Y" || $utr == "y" || $utr == "Yes" ]]; then
        print_fluorescent_yellow "    11.  Adding UTRs..."
        python3 $pmetroot/parse_utrs.py      \
            $indexingOutputDir/promoters.bed \
            $indexingOutputDir/sorted.gff3   \
            $universefile
    else
        print_fluorescent_yellow "    11.  (skipped) Adding UTRs..."
    fi

    # -------------------------------------------------------------------------------------------
    # 12. promoter lenfths from promoters.bed
    print_fluorescent_yellow "    12.  Promoter lengths from promoters.bed (promoter_lengths_all.txt)"
    # python3 $pmetroot/parse_promoter_lengths.py \
    #     $indexingOutputDir/promoters.bed \
    #     $indexingOutputDir/promoter_lengths.txt
    awk '{print $4 "\t" ($3 - $2)}'      \
        $indexingOutputDir/promoters.bed \
        > $indexingOutputDir/promoter_lengths.txt

    # -------------------------------------------------------------------------------------------
    # 13. Update genes list
    print_fluorescent_yellow "    13.  Update genes list: complete list of all genes found (universe.txt)"
    cut -f 1 $indexingOutputDir/promoter_lengths.txt > $universefile

    # -------------------------------------------------------------------------------------------
    # 14. create promoters fasta
    print_fluorescent_yellow "    14.  Creating promoters file (promoters_rough.fa)";
    bedtools getfasta -fi                         \
        $indexingOutputDir/genome_stripped.fa     \
        -bed $indexingOutputDir/promoters.bed     \
        -fo $indexingOutputDir/promoters_rough.fa \
        -name

    # -------------------------------------------------------------------------------------------
    # 15. replace the id of each seq with gene names
    print_fluorescent_yellow "    15.  Replacing the id of each sequences' with gene names (promoters.fa)"
    sed 's/::.*//g' $indexingOutputDir/promoters_rough.fa > $indexingOutputDir/promoters.fa

    # -------------------------------------------------------------------------------------------
    # 16. promoters.bg from promoters.fa
    print_fluorescent_yellow "    16.  fasta-get-markov: a Markov model from promoters.fa. (promoters.bg)"
    fasta-get-markov $indexingOutputDir/promoters.fa > $indexingOutputDir/promoters.bg

    # -------------------------------------------------------------------------------------------
    # 17. individual motif files from user's meme file
    print_fluorescent_yellow "    17.  Spliting motifs into individual meme files (folder memefiles)"
    mkdir -p $indexingOutputDir/memefiles
    python3 $pmetroot/parse_memefile_batches.py $memefile $indexingOutputDir/memefiles/ $threads

    # -------------------------------------------------------------------------------------------
    # 18. IC.txt
    print_fluorescent_yellow "    18.  Generating information content (IC.txt)"
    mkdir -p $indexingOutputDir/memefilestemp
    python3 $pmetroot/parse_memefile.py $memefile $indexingOutputDir/memefilestemp/
    python3 $pmetroot/calculateICfrommeme_IC_to_csv.py \
        $indexingOutputDir/memefilestemp/              \
        $indexingOutputDir/IC.txt
    rm -rf $indexingOutputDir/memefilestemp/

    if [ true ]; then
        # rm -rf $indexingOutputDir/IC.txt
        rm -rf $indexingOutputDir/bedgenome.genome
        # rm -rf $indexingOutputDir/fimohits
        rm -rf $indexingOutputDir/genelines.bed
        rm -rf $indexingOutputDir/genelines.gff3
        rm -rf $indexingOutputDir/genome_stripped.fa
        rm -rf $indexingOutputDir/genome_stripped.fa.fai
        # rm -rf $indexingOutputDir/histogram_distance_tss.png
        # rm -rf $indexingOutputDir/invalid_genelines.txt
        rm -rf $indexingOutputDir/invalid_gff3_lines.txt
        # rm -rf $indexingOutputDir/length_to_tss.txt
        # rm -rf $indexingOutputDir/memefiles
        # rm -rf $indexingOutputDir/promoter_lengths.txt
        rm -rf $indexingOutputDir/promoters.bed
        rm -rf $indexingOutputDir/promoters_rough.fa
        rm -rf $indexingOutputDir/sorted.gff3
    fi

    for maxk in ${maxkRange[@]}; do
        for topn in ${topnRange[@]}; do
            cp -r $indexingOutputDir $outputDir/TEMP_LEN${promlength}_K${maxk}_N${topn}_FIMO${fimothresh//./}
        done
    done

    rm -rf $indexingOutputDir
done

# 计算总迭代次数 total number of iterations
totalIterations=$(( ${#promlengthRange[@]} * ${#maxkRange[@]} * ${#topnRange[@]} ))
# 初始化已完成迭代次数的计数器 Initialization of the iteration counter has been completed.
completedIterations=0

for promlength in ${promlengthRange[@]}; do
    for maxk in ${maxkRange[@]}; do
        for topn in ${topnRange[@]}; do
            indexingOutputDir=$outputDir/TEMP_LEN${promlength}_K${maxk}_N${topn}_FIMO${fimothresh//./}

            # 更新已完成迭代次数 Update the number of completed iterations
            ((completedIterations++))
            # -------------------------------- Run fimo and pmetindex --------------------------
            print_green "$(date '+%Y-%m-%d %H:%M:%S') - (${completedIterations}/${totalIterations}) Running FIMO and PMET indexing with ${threads} CPUS......"
            print_green "    Promoter legnth: $promlength"
            print_green "    MaxK           : $maxk"
            print_green "    Topn           : $topn"
            print_green "    Fimo threshold : $fimothresh"
            print_green "    Indexing output: $indexingOutputDir"

            start_time=$SECONDS

            runFimoIndexing () {
                memefile=$1
                indexingOutputDir=$2
                fimothresh=$3
                pmetroot=$4
                maxk=$5
                topn=$6

                $pmetroot/fimo                               \
                    --topk $maxk                             \
                    --topn $topn                             \
                    --text                                   \
                    --no-qvalue                              \
                    --thresh $fimothresh                     \
                    --verbosity 1                            \
                    --oc $indexingOutputDir/fimohits         \
                    --bgfile $indexingOutputDir/promoters.bg \
                    $memefile                                \
                    $indexingOutputDir/promoters.fa          \
                    $indexingOutputDir/promoter_lengths.txt
            }
            export -f runFimoIndexing

            numfiles=$(ls -l $indexingOutputDir/memefiles/*.txt | wc -l)
            print_orange "    $numfiles motif meme files found"

            find $indexingOutputDir/memefiles -name \*.txt \
                | parallel --bar --jobs=$threads \
                    "runFimoIndexing {} $indexingOutputDir $fimothresh $pmetroot $maxk $topn"

            mv $indexingOutputDir/fimohits/binomial_thresholds.txt $indexingOutputDir/

            # when indexing completed, remove TEMP flag from folder
            mv $indexingOutputDir $outputDir/LEN${promlength}_K${maxk}_N${topn}_FIMO${fimothresh//./}


            end_time=$SECONDS
            elapsed_time=$((end_time - start_time))
            days=$((elapsed_time/86400))
            hours=$(( (elapsed_time%86400)/3600 ))
            minutes=$(( (elapsed_time%3600)/60 ))
            seconds=$((elapsed_time%60))
            print_orange "(${completedIterations}/${totalIterations}) FIMO and PMET indexing taken: $days day $hours hour $minutes minute $seconds second\n"
        done
    done
done

print_green "DONE"
