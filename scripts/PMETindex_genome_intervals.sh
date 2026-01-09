#!/bin/bash
set -e

# 22.1.18 Charlotte Rich
# last edit: 7.2.18 - removed the make 1 big fimohits files

# cl_index_wrapper.sh
# mac -> server Version differences
# ggrep = grep

#Called when user selects 'Genomic Intervals'
#Input files are genomic intevals fasta file, meme file location, gene clusters file
#Other inputs N and k


function usage () {
    cat >&2 <<EOF
USAGE: PMET_index_ATAC_fasta_peaks.sh [options] <genomefile> <peaks.bed> <memefile>

Creates PMET index for Paired Motif Enrichment Test from ATAC seq data using peaks in fasta and bed files.
Required arguments:
-r <PMETindex_path>	: Full path of PMET_index. Required.

Optional arguments:
-o <output_directory> : Output directory for results
-n <topn>	: How many top promoter hits to take per motif. Default=5000
-k <max_k>	: Maximum motif hits allowed within each promoter.  Default: 5
-j <max_jobs>	: Max number of jobs to run at once.  Default: 1

EOF
}

function error_exit() {
  echo "ERROR: $1" >&2
  usage
  exit 1
}


# set up arguments
topn=5000
maxk=5
fimothresh=0.05
pmetroot="scripts"
progFile="progress/progress"

outputdir=
genomefile=
memefile=
gene_input_file=

# check if arguments have been specified
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"  >&2
    usage
    exit 1
fi

# bring in arguments
while getopts ":r:o:k:n:g:f" opt; do
  case $opt in
    r) echo "Full path of PMET_index:  $OPTARG" >&2
    pmetroot=$OPTARG;;
    o) echo "Output directory for results: $OPTARG" >&2
    outputdir=$OPTARG;;
    n) echo "Top n promoter hits to take per motif: $OPTARG" >&2
    topn=$OPTARG;;
    k) echo "Top k motif hits within each promoter: $OPTARG" >&2
    maxk=$OPTARG;;
    f) echo "Fimo threshold: $OPTARG" >&2
    fimothresh=$OPTARG;;
    g) echo "Progress file: $OPTARG" >&2
    progFile=$OPTARG;;
    \?) echo "Invalid option: -$OPTARG" >&2
    exit 1;;
    :)  echo "Option -$OPTARG requires an argument." >&2
    exit 1;;
  esac
done

#rename input file variable
shift $((OPTIND - 1))


genomefile=$1
memefile=$2
gene_input_file=$3

[ ! -d $outputdir ] && mkdir $outputdir
cd $outputdir

start=$SECONDS

# progress taken from 1% to 95%
echo -e "0.01\tPreparing sequences..." > $progFile
echo "Preparing sequences...";




#final pmet binary requires the universe file. Need to create this if validation scrip didnt.
#In promoters version, this is initially all genes in gff3 file. This version is used to add UTRs if
#requested, but any genes not in promoter_lengths file are filtered out before we get to PMET binary stage
#In this version we can just take a copy of all IDs in promoter lengths as we dont to UTR stuff


universefile=$outputdir/universe.txt

if [[ ! -f "$universefile" || ! -f "$outputdir/promoter_lengths.txt" ]]; then

	# should have been done by consistency checker
	# *** ADD THE DEPUPLICATION OF THE FASTA FILE HERE ****
	python3 $pmetroot/python/deduplicate.py $genomefile $outputdir/no_duplicates.fa

	# generate the promoter lengths file from the fasta file
	python3 $pmetroot/python/parse_promoter_lengths_from_fasta.py $outputdir/no_duplicates.fa $outputdir/promoter_lengths.txt
	rm -f $outputdir/no_duplicates.fa

  cut -f 1  $outputdir/promoter_lengths.txt > $universefile
fi



#now we can actually FIMO our way to victory
/usr/local/meme/bin/fasta-get-markov $genomefile > $outputdir/genome.bg
#FIMO barfs ALL the output. that's not good. time for individual FIMOs
#on individual MEME-friendly motif files too

duration=$(( SECONDS - start ))
echo $duration" secs"

start=$SECONDS

echo "Processing motifs...";
echo -e "0.05\tProcessing motifs..." >$progFile

###†Make motif  files from user's meme file
[ ! -d $outputdir/memefiles ] && mkdir $outputdir/memefiles

python3 $pmetroot/python/parse_memefile.py $memefile $outputdir/memefiles/

### creates IC.txt tsv file from, motif files
python3 $pmetroot/python/calculateICfrommeme_IC_to_csv.py $outputdir/memefiles/ $outputdir/IC.txt

###†Create a fimo hits file form each motif file

[ ! -d $outputdir/fimo ] && mkdir $outputdir/fimo
[ ! -d $outputdir/fimohits ] && mkdir $outputdir/fimohits

shopt -s nullglob #prevent loop produncing '*.txt'
#20 sec min per file, progress goes from 5% to 70%

numfiles=$(ls -l $outputdir/memefiles/*.txt | wc -l)
echo $numfiles" found"
n=0
#paralellise this loop
for memefile in $outputdir/memefiles/*.txt; do
    let n=$n+1
    prog="0"$(bc <<< "scale=2;0.05+0.65*$n/$numfiles")
    echo -e "$prog\tProcessing motif $n of $numfiles" > $progFile
    echo "$prog\tProcessing motif $n of $numfiles"
    sleep 0.1
    bfid=`basename $memefile`
    fimofile=$bfid

    /usr/local/meme/bin/fimo --text --thresh $fimothresh --verbosity 1 --bgfile $outputdir/genome.bg $memefile $genomefile > $outputdir/fimo/$fimofile

done

################
rm -r $outputdir/memefiles
rm $outputdir/genome.bg

duration=$(( SECONDS - start ))
echo $duration" secs"

start=$SECONDS

#final run indexing. This needs N, k , location of fimo files and the output directory

#sorts and filters each fimo/XXX.txt file according to k and N and writes output to fimohits/XXX.txt
#created binomial thresholds as it does this
#8 sec per file, takes progress from 70% to 95%
echo -e "0.7\tIndexing..." > $progFile
echo "Indexing...";
sleep 0.1

#promoter_lengths file must have an entry for every gene in fimo files

$pmetroot/pmetindex -f $outputdir/fimo -k $maxk -n $topn -o $outputdir -p $outputdir/promoter_lengths.txt -g $progFile
#equivalent of parse_matrix_n.py
#This creates binomial_thresholds.txt
echo -e "0.95\tDone" > $progFile
echo "Done";
duration=$(( SECONDS - start ))
echo $duration" secs"


# Match lines from genes file to any string in pattern, which is content of universe file
grep -Ff $universefile $gene_input_file > $gene_input_file"temp"
mv $gene_input_file"temp" $gene_input_file


exit 0;


#next stage needs the following inputs

#   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
#   bimnomial_thresholds.txt    made by PMETindex
#   IC.txt                      made by calculateICfrommeme.py from meme file
#   gene input file             supplied by user
