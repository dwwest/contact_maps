import numpy as np

lines = []
with open('GSM2418860_WT_CTCF_peaks.txt', 'r') as f:
    for ind, line in enumerate(f):
        if ind != 0:
            lines.append('\t'.join(line.strip('\n').split('\t')[0:3])+'\n')

with open('GSM2418860_WT_CTCF_peaks.bed', 'w+') as f:
    for line in lines:
        f.write(line)

