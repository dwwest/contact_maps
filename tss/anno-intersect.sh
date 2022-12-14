#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate microc

sample=$1
bed=$2

bam="${sample}.PT.bam"
pair="$sample-filt.pairs"
anno=`basename $bed .bed`

# make new bams with only reads that fall in heterochromatin or euchromatin regions 
bedtools intersect -a $bam -b $bed -wa | samtools view > $anno-temp.sam

# filter .pairs files to include only lines from anno bam (matched by names) to make
#       anno .pairs file for contact probability analysis
awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' $anno-temp.sam $pair > $anno-$pair

rm "$anno-temp.sam"

