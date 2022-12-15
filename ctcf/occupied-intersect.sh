#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate microc

bedtools intersect -wa -a granges_expanded.bed -b GSM2418860_WT_CTCF_peaks.bed > intersected_GSM2418860_expanded.bed
