#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate bayes

bigWigToWig GSM2418860_WT_CTCF.bw GSM2418860_WT_CTCF.wig

wig2bed -x < GSM2418860_WT_CTCF.wig > GSM2418860_WT_CTCF.bed

