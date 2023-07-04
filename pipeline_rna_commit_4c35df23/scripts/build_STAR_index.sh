#!/bin/bash

usage="Usage: $0 -S sjdbOverhang -O output_directory [-R reference -A annotation -C nof_cores]"

base_dir="/path/to/references/Hsap_GRCh37/RNA/"

nof_cores=1
reference="$base_dir/genome.hs37d5.fa"
annotation="$base_dir/gencode.v19.annotation.hs37d5_chr.gtf"

while getopts R:A:S:O:C: flag; do
    case $flag in
    R)
    reference=$OPTARG
    ;;
    A)
    annotation=$OPTARG
    ;;
    S)
    sjdbOverhang=$OPTARG
    ;;
    O)
    output_dir=$OPTARG
    ;;
    C)
    nof_cores=$OPTARG
    ;;
    \?)
    echo $usage >& 2
    exit -1
    ;;
    esac
done
shift $(( OPTIND - 1));


if [[ -z $sjdbOverhang ]]; then
    echo "Need to pass sjdbOverhang to ${0}. Exiting ... "
    echo $usage
    exit -1
fi


if [[ -z $output_dir ]]; then
    echo "Need to pass output directory to ${0}. Exiting ... "
    echo $usage
    exit -1
fi

   
if [[ ! -d $output_dir ]]; then
    mkdir -p $output_dir
fi  


master_time_start=`date +%s`
echo ">>>>>>>>>" >&2
echo Start: `hostname --short` `date` >&2  

echo "Running STAR genomeGenerate using the following parameters:" >&2
echo "--genomeDir $output_dir" >&2
echo "--genomeFastaFiles $reference" >&2
echo "--sjdbOverhang $sjdbOverhang" >&2
echo "--sjdbGTFfile $annotation" >&2
echo "--runThreadN $nof_cores" >&2

  
STAR \
--runMode genomeGenerate \
--genomeDir $output_dir \
--genomeFastaFiles $reference \
--sjdbOverhang $sjdbOverhang \
--sjdbGTFfile $annotation \
--runThreadN $nof_cores



master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))` 
echo "Finished running STAR genomeGenerate" >&2
echo "STAR genomeGenerate completed in $master_time_exec seconds") >&2
echo End: `hostname --short` `date` >&2
