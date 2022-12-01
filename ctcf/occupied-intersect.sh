#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate microc

bedtools intersect -a GSM2418860_WT_CTCF_peaks.bed -b granges_expanded.bed > intersected_GSM2418860_expanded.bed
