#!/bin/bash
set -e

# 22.1.18 Charlotte Rich
# last edit: 7.2.18 - removed the make 1 big fimohits files

# cl_index_wrapper.sh
# mac -> server Version differences
# ggrep = grep


function usage () {
    cat >&2 <<EOF
USAGE: PMET_index_ATAC_fasta_peaks.sh [options] <peaks.fa> <peaks.bed> <memefile>

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

function check_file() {
  file=$1
  name=$2
  flag=$3

  if [ ! -f "$file" ]
  then error_exit "$name is not a file.  $flag flag is required"
fi
}

function check_dir() {
  dir=$1
  name=$2
  flag=$3

  if [ ! -d "$dir" ]
  then error_exit "$name is not a directory.  $flag flag is required"
fi
}


function check_set() {
  value=$1
  name=$2
  flag=$3

  # -z string  = True if the length of string is zero.
  # exits if a string has not been specified
  if [ -z "$value" ]
  then error_exit "$name has not been specified.  $flag flag is required"
fi

}
# set up arguments
topn=5000
maxk=5
maxjobs=1
# set up empty variables
pmetroot=
outputdir=
peaksfasta=
peaksbed=
memefile=

# check if arguments have been specified
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"  >&2
    usage
    exit 1
fi

# bring in arguments
while getopts ":r:o:k:n:j:" opt; do
  case $opt in
    r) echo "Full path of PMET_index:  $OPTARG" >&2
    pmetroot=$OPTARG;;
    o) echo "Output directory for results: $OPTARG" >&2
    outputdir=$OPTARG;;
    n) echo "Top n promoter hits to take per motif: $OPTARG" >&2
    topn=$OPTARG;;
    k) echo "Top k motif hits within each promoter: $OPTARG" >&2
    maxk=$OPTARG;;
    j) echo "Max no. parallel jobs: $OPTARG" >&2
    maxjobs=$OPTARG;;
    \?) echo "Invalid option: -$OPTARG" >&2
    exit 1;;
    :)  echo "Option -$OPTARG requires an argument." >&2
    exit 1;;
  esac
done

#rename input file variable
shift $((OPTIND - 1))
date

peaksfasta=$1
peaksbed=$2
memedir=$3

# check things exist
check_set "$pmetroot" "PMET index scripts path" "-r"
check_file "$peaksfasta" "peak fasta input file" ""
check_file "$peaksbed" "peaks bed file" ""
check_file "$memedir" "meme input file" ""

# check if output directory is specified, if it is check it already exists, if not create it
# and enter that directory
# if a motif_found_final.txt file is already in the folder exit to avoid overwriting data.
if [ -n $outputdir ]; then
echo "Output directory:" $outputdir
  if [ -d $outputdir ]; then
    cd $outputdir
    if [ -d fimo ]; then
      error_exit "Output directory already exists and contains PMET_index data, exiting to avoid overwriting"
    fi
  else
    echo "Output directory does not exist, creating directory"
    mkdir $outputdir
    cp $peaksfasta $outputdir
    cp $peaksbed $outputdir
    cp -r $memedir $outputdir
    cd $outputdir
  fi
fi


# samtools faidx genome_stripped.fa
if [ ! -f peaks.fa ]
then
  cp $peaksfasta peaks.fa
fi

if [ ! -f peaks.bed ]
then
  cp $peaksbed peaks.bed
fi

# get peak lengths
#python3 $pmetroot/parse_promoter_lengths.py peaks.bed
python3 $pmetroot/python/parse_promoter_lengths_from_fasta.py peaks.fa

#now we can actually FIMO our way to victory
fasta-get-markov peaks.fa > peaks.bg
#FIMO barfs ALL the output. that's not good. time for individual FIMOs
#on individual MEME-friendly motif files too

# this just copies the directory rather than splitting up one memefile, leave for now
#if [ ! -d memefiles ]
#then
#  cp -r $memedir memefiles
#fi

mkdir memefiles

python3 $pmetroot/python/parse_memefile.py $memedir

python3 $pmetroot/python/calculateICfrommeme.py


runIndexing () {
  fid=$1
  echo $fid
  bfid=`basename $fid`
  echo $bfid
  progpath=$4
  # get all the possible motif hits using fimo
  fimo --text --thresh 0.05 --verbosity 1 --bgfile peaks.bg $fid peaks.fa > fimo_$bfid
  # # parse the fimo output to get top n promoters containing top n hits
  python3 $progpath/parse_matrix_n.py fimo_$bfid $2 $3
  rm fimo_$bfid
}
export -f runIndexing

mkdir fimohits
find memefiles -name \*.txt | parallel --jobs=$maxjobs "runIndexing {} $maxk $topn $pmetroot"
wait



#there's a lot of intermediate files that need blanking

rm fimo_*
#rm peaks.bed
rm peaks.bg
rm peaks.fa
date
