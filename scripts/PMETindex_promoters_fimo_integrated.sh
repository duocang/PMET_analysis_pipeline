#!/bin/bash
set -e

# Load shared utilities
script_dir=$(cd -- "$(dirname "$0")" && pwd)
source "$script_dir/lib/print_colors.sh"
source "$script_dir/lib/timer.sh"

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

# Override error_exit to show usage
error_exit() { echo "ERROR: $1" >&2; usage; exit 1; }


# set up defaults
topn=5000
maxk=5
promlength=1000
fimothresh=0.05
overlap="AllowOverlap"
utr="No"
gff3id='gene_id'
toolDir="scripts"
pmetDir="build"
threads=4
delete=yes
isPoisson=no

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

while getopts ":r:i:o:n:k:p:f:v:u:t:d:x:" options; do
    case $options in
        r) toolDir=$OPTARG;;
        i) gff3id=$OPTARG;;
        o) indexingOutputDir=$OPTARG;;
        n) topn=$OPTARG;;
        k) maxk=$OPTARG;;
        p) promlength=$OPTARG;;
        f) fimothresh=$OPTARG;;
        v) overlap=$OPTARG;;
        u) utr=$OPTARG;;
        t) threads=$OPTARG;;
        d) delete=$OPTARG;;
        x) isPoisson=$OPTARG;;
        \?) print_red  "Invalid option: -$OPTARG" >&2
        exit 1;;
        :)  print_red "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done

shift $((OPTIND - 1))
genomefile="$1"
gff3file="$2"
memefile="$3"
universefile="$indexingOutputDir/universe.txt"
bedfile="$indexingOutputDir/genelines.bed"

print_white "Genome file                  : "; print_orange "$genomefile"
print_white "Annotation file              : "; print_orange "$gff3file"
print_white "Motif meme file              : "; print_orange "$memefile"

print_white "PMET index path              : "; print_orange "$toolDir"
print_white "GFF3 identifier              : "; print_orange "$gff3id"
print_white "Output directory             : "; print_orange "$indexingOutputDir"
print_white "Top n promoters              : "; print_orange "$topn"  # Default to 5000 if not set
print_white "Top k motif hits             : "; print_orange "$maxk"     # Default to 5 if not set
print_white "Length of promoter           : "; print_orange "$promlength"  # Default to 1000 if not set
print_white "Fimo threshold               : "; print_orange "$fimothresh"
print_white "Promoter overlap handling    : "; print_orange "$overlap"
print_white "Include 5' UTR               : "; print_orange "$utr"
print_white "Number of threads            : "; print_orange "$threads"

mkdir -p $indexingOutputDir

start_time=$SECONDS
print_green "Preparing data for FIMO and PMET index..."
# -------------------------------------------------------------------------------------------
# 1. sort annotaion by gene coordinates
print_fluorescent_yellow "     1. Sorting annotation by gene coordinates"
chmod a+x "$toolDir/gff3sort/gff3sort.pl"
"$toolDir/gff3sort/gff3sort.pl" "$gff3file" > "$indexingOutputDir/sorted.gff3"

# -------------------------------------------------------------------------------------------
# 2. extract gene line from annoitation
print_fluorescent_yellow "     2. Extracting gene line from annoitation"
if [[ "$(uname)" == "Linux" ]]; then
    grep -P '\tgene\t' "$indexingOutputDir/sorted.gff3" > "$indexingOutputDir/genelines.gff3"
elif [[ "$(uname)" == "Darwin" ]]; then
    grep '\tgene\t' "$indexingOutputDir/sorted.gff3" > "$indexingOutputDir/genelines.gff3"
else
    print_red "Unsupported operating system."
    exit 1
fi

# -------------------------------------------------------------------------------------------
# 3. extract chromosome , start, end, gene ('gene_id' for input)
print_fluorescent_yellow "     3. Extracting chromosome, start, end, gene"
# parse up the .bed for promoter extraction, 'gene_id'
# 使用grep查找字符串 check if gene_id is present
if grep -q "$gff3id" "$indexingOutputDir/genelines.gff3"; then
    python3 "$toolDir/parse_genelines.py" "$gff3id" "$indexingOutputDir/genelines.gff3" "$bedfile"
else
    gff3id='ID='
    python3 "$toolDir/parse_genelines.py" "$gff3id" "$indexingOutputDir/genelines.gff3" "$bedfile"
fi
# -------------------------------------------------------------------------------------------
# 4. filter invalid genes: start should be smaller than end
print_fluorescent_yellow "     4. Filter invalid coordinates: start > end"
awk '$2 >= $3' "$bedfile" > "$indexingOutputDir/invalid_gff3_lines.txt"
awk '$2 <  $3' "$bedfile" > "$indexingOutputDir/genelines_valid.bed"
mv "$indexingOutputDir/genelines_valid.bed" "$bedfile"
# LC_ALL=C sort -k1,1 -k6,6 -k2,2n "$bedfile" -o "$bedfile"
# 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
# In the BED file format, the start position is always less than the end position for both positive (+) and negative (-) chains.
# 起始和终止位置是指定基因上的物理位置，而不是表达或翻译的方向。
# start and end positions specify the physical location of the gene, rather than the direction of expression or translation.

print_fluorescent_yellow "        Calculate length of chromosome (length_of_chromosome.txt)"
python3 "$toolDir/calculate_chromosome_length.py" "$gff3file" "$indexingOutputDir/length_of_chromosome.txt"

print_fluorescent_yellow "        Calculate length of space to TSS (length_to_tss.txt)"
python3 "$toolDir/calculate_length_to_tss.py" \
    "$bedfile" \
    "$indexingOutputDir/length_of_chromosome.txt" \
    "$indexingOutputDir/length_to_tss.txt"

# -------------------------------------------------------------------------------------------
# 5. list of all genes found
print_fluorescent_yellow "\n     5. Extracting genes names: complete list of all genes found (universe.txt)"
cut -f 4 "$bedfile" > "$universefile"

# -------------------------------------------------------------------------------------------
# 6. strip the potential FASTA line breaks. creates genome_stripped.fa
print_fluorescent_yellow "     6. Removing potential FASTA line breaks (genome_stripped.fa)"
awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
    { printf "%s",$0;} \
    END { print ""; }' "$genomefile" > "$indexingOutputDir/genome_stripped.fa"

# -------------------------------------------------------------------------------------------
# 7. create the .genome file which contains coordinates for each chromosome start
print_fluorescent_yellow "     7. Listing chromosome start coordinates (bedgenome.genome)"
samtools faidx "$indexingOutputDir/genome_stripped.fa"
cut -f 1-2 "$indexingOutputDir/genome_stripped.fa.fai" > "$indexingOutputDir/bedgenome.genome"

# -------------------------------------------------------------------------------------------
# 8. create promoters' coordinates from annotation
print_fluorescent_yellow "     8. Creating promoters' coordinates from annotation (promoters.bed)"
# 在bedtools中，flank是一个命令行工具，用于在BED格式的基因组坐标文件中对每个区域进行扩展或缩短。
# In bedtools, flank is a command-line tool used to extend or shorten each region in a BED format genomic coordinate file.
# 当遇到负链（negative strand）时，在区域的右侧进行扩展或缩短，而不是左侧。
# When a negative strand is encountered, it is expanded or shortened on the right side of the region, not the left.
bedtools flank                             \
    -l "$promlength"                       \
    -r 0 -s -i "$bedfile"                  \
    -g "$indexingOutputDir/bedgenome.genome" \
    > "$indexingOutputDir/promoters_not_sorted.bed"
# Sort by starting coordinate
sortBed -i "$indexingOutputDir/promoters_not_sorted.bed" > "$indexingOutputDir/promoters.bed"
rm -rf "$indexingOutputDir/promoters_not_sorted.bed"

# -------------------------------------------------------------------------------------------
# 9. remove overlapping promoter chunks
overlap_lower=$(echo "$overlap" | tr '[:upper:]' '[:lower:]')
if [[ "$overlap_lower" == "nooverlap" || "$overlap_lower" == "no" || "$overlap_lower" == "n" ]]; then
	print_fluorescent_yellow "     9. Removing overlapping promoter chunks (promoters.bed)"
	sleep 0.1
	bedtools subtract                       \
        -a "$indexingOutputDir/promoters.bed" \
        -b "$bedfile"                         \
        > "$indexingOutputDir/promoters2.bed"
	mv "$indexingOutputDir/promoters2.bed" "$indexingOutputDir/promoters.bed"
else
    print_fluorescent_yellow "     9. (skipped) Removing overlapping promoter chunks (promoters.bed)"
fi


# -------------------------------------------------------------------------------------------
# 10. check split promoters. if so, keep the bit closer to the TSS
print_fluorescent_yellow "    10. Checking split promoter (if so):  keep the bit closer to the TSS (promoters.bed)"
python3 "$toolDir/assess_integrity.py" "$indexingOutputDir/promoters.bed"

# -------------------------------------------------------------------------------------------
# 11. add 5' UTR
utr_lower=$(echo "$utr" | tr '[:upper:]' '[:lower:]')
if [[ "$utr_lower" == "yes" || "$utr_lower" == "y" ]]; then
    print_fluorescent_yellow "    11. Adding UTRs...";
	python3 "$toolDir/parse_utrs.py"      \
        "$indexingOutputDir/promoters.bed" \
        "$indexingOutputDir/sorted.gff3"   \
        "$universefile"
else
    print_fluorescent_yellow "    11. (skipped) Adding UTRs...";
fi

# -------------------------------------------------------------------------------------------
# 12. promoter lengths from promoters.bed
print_fluorescent_yellow "    12. Promoter lengths from promoters.bed (promoter_lengths.txt)"
awk '{print $4 "\t" ($3 - $2)}' "$indexingOutputDir/promoters.bed" \
    > "$indexingOutputDir/promoter_lengths.txt"

# -------------------------------------------------------------------------------------------
# 13. Update genes list
print_fluorescent_yellow "    13. Update genes list: complete list of all genes found (universe.txt)"
cut -f 1 "$indexingOutputDir/promoter_lengths.txt" > "$universefile"

# -------------------------------------------------------------------------------------------
# 14. create promoters fasta
print_fluorescent_yellow "    14. Creating promoters file (promoters_rough.fa)";
bedtools getfasta -fi \
    "$indexingOutputDir/genome_stripped.fa"     \
    -bed "$indexingOutputDir/promoters.bed"     \
    -fo "$indexingOutputDir/promoters_rough.fa" \
    -name

# -------------------------------------------------------------------------------------------
# 15. replace the id of each seq with gene names
print_fluorescent_yellow "    15. Replacing id of each sequences' with gene names (promoters.fa)"
sed 's/::.*//g' "$indexingOutputDir/promoters_rough.fa" > "$indexingOutputDir/promoters.fa"

# -------------------------------------------------------------------------------------------
# 16. promoters.bg from promoters.fa
print_fluorescent_yellow "    16. fasta-get-markov: a Markov model from promoters.fa. (promoters.bg)"
fasta-get-markov "$indexingOutputDir/promoters.fa" > "$indexingOutputDir/promoters.bg"

# -------------------------------------------------------------------------------------------
# 17. IC.txt
print_fluorescent_yellow "    17. Generating information content (IC.txt)"
mkdir -p "$indexingOutputDir/memefiles"
python3 "$toolDir/parse_memefile.py" "$memefile" "$indexingOutputDir/memefiles/"
python3                                         \
    "$toolDir/calculateICfrommeme_IC_to_csv.py" \
    "$indexingOutputDir/memefiles/"             \
    "$indexingOutputDir/IC.txt"

# -------------------------------------------------------------------------------------------
# 18. individual motif files from user's meme file (reusing memefiles directory from step 17)
print_fluorescent_yellow "    18. Spliting motifs into individual meme files (folder memefiles)"
rm -rf "$indexingOutputDir/memefiles/"*
python3 "$toolDir/parse_memefile_batches.py" "$memefile" "$indexingOutputDir/memefiles/" "$threads"

# -------------------------------- Run fimo and pmetindex --------------------------
# mkdir -p $indexingOutputDir/fimohits

print_green "Running FIMO and PMET index..."
runFimoIndexing () {
    local memefile="$1"
    local indexingOutputDir="$2"
    local fimothresh="$3"
    local pmetDir="$4"
    local maxk="$5"
    local topn="$6"
    local isPoisson="$7"
    local isPoisson_lower=$(echo "$isPoisson" | tr '[:upper:]' '[:lower:]')
    local poisson_flag=""

    if [[ "$isPoisson_lower" == "true" || "$isPoisson_lower" == "t" || "$isPoisson_lower" == "yes" || "$isPoisson_lower" == "y" ]]; then
        poisson_flag="--poisson"
    fi

    "$pmetDir/index_fimo_fused"                     \
        $poisson_flag                               \
        --no-qvalue                                 \
        --text                                      \
        --thresh "$fimothresh"                      \
        --verbosity 1                               \
        --bgfile "$indexingOutputDir/promoters.bg"  \
        --topn "$topn"                              \
        --topk "$maxk"                              \
        --oc "$indexingOutputDir"                   \
        "$memefile"                                 \
        "$indexingOutputDir/promoters.fa"           \
        "$indexingOutputDir/promoter_lengths.txt"
}
export -f runFimoIndexing

nummotifs=$(grep -c '^MOTIF' "$memefile")
print_orange "    $nummotifs motifs found"

find "$indexingOutputDir/memefiles" -name "*.txt" \
    | parallel --bar --jobs="$threads" \
        runFimoIndexing {} "$indexingOutputDir" "$fimothresh" "$pmetDir" "$maxk" "$topn" "$isPoisson"

# mv "$indexingOutputDir/fimohits/binomial_thresholds.txt" "$indexingOutputDir/"

delete_lower=$(echo "$delete" | tr '[:upper:]' '[:lower:]')
if [[ "$delete_lower" == "yes" || "$delete_lower" == "y" ]]; then
    print_green "Deleting unnecessary files..."
    rm -rf \
        "$indexingOutputDir/bedgenome.genome" \
        "$indexingOutputDir/genelines.bed" \
        "$indexingOutputDir/genelines.gff3" \
        "$indexingOutputDir/genome_stripped.fa" \
        "$indexingOutputDir/genome_stripped.fa.fai" \
        "$indexingOutputDir/memefiles" \
        "$indexingOutputDir/promoters.bed" \
        "$indexingOutputDir/promoters.bg" \
        "$indexingOutputDir/promoters.fa" \
        "$indexingOutputDir/promoters_rough.fa" \
        "$indexingOutputDir/sorted.gff3"
fi

# 计算 $indexingOutputDir/fimohits 目录下 .txt 文件的数量
# Count the number of .txt files in the $indexingOutputDir/fimohits directory
file_count=$(find "$indexingOutputDir/fimohits" -maxdepth 1 -type f -name "*.txt" | wc -l)

# 检查文件数量是否等于 motif 的数量 （$nummotifs）
# Check if the number of files equals the number of motifs ($nummotifs)
if [ "$file_count" -eq "$nummotifs" ]; then
    print_elapsed_time $start_time
else
    print_elapsed_time_with_status $start_time "error" "there are $file_count fimohits files, it should be $nummotifs."
fi

# next stage needs the following inputs

#   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
#   bimnomial_thresholds.txt    made by PMETindex
#   IC.txt                      made by calculateICfrommeme.py from meme file
#   gene input file             supplied by user
