import numpy as np
import matplotlib.pyplot as plt
import sys
import os

"""
Maps pairs to ctcf sites for making a contact map
"""

# this should be twice + 20 whatever the value is in expand_granges.py
chromatin_len = 1600

slice_one = sys.argv[1]
slice_two = sys.argv[2]
last = sys.argv[3]

pairs = []
with open("granges_expanded-pubWTmESC-mapped-filt.pairs") as f:
    for line in f:
        current = line.strip('\n').split('\t')
        # filter out any pairs that are more than chromatin_len bp apart
        # because we only want pairs from our ROI
        if abs(float(current[2]) - float(current[4])) <= chromatin_len:
            pairs.append([])
            pairs[-1].append([current[1], current[5], int(current[2])])
            pairs[-1].append([current[3], current[6], int(current[4])])
            pairs[-1].append(current[7])
if last == False:
    pairs = pairs[slice_one:slice_two]
elif last == True:
    pairs = pairs[slice_one:]
print('Pairs loaded...')

ctcf = {}
with open("granges_expanded.bed") as f:
    for ind,line in enumerate(f):
        if ind > 0:
            current = line.strip('\n').split('\t')
            if current[0] not in ctcf.keys():
                ctcf[current[0]] = [int(current[1])]
            else:
                ctcf[current[0]].append(int(current[1]))
print('Sites loaded...')

current_chrom = ''
contact_map = np.zeros([chromatin_len+1, chromatin_len+1, 4])
for pair_ind, pair in enumerate(pairs):
    chrom_one = pairs[pair_ind][0][0]
    site_one = pairs[pair_ind][0][2]
    strand_one = pairs[pair_ind][0][1]
    chrom_two = pairs[pair_ind][1][0]
    site_two = pairs[pair_ind][1][2]
    strand_two = pairs[pair_ind][1][1]
    if chrom_one == chrom_two:
        if chrom_one != current_chrom:
            current_chrom = chrom_one
            print(current_chrom)
        for s in ctcf[chrom_one]:
            if site_one >= s and site_one <= chromatin_len + s:
                if site_two >= s and site_two <= chromatin_len + s:
                    if abs(site_one - site_two) >= 147:
                        """ INWARD """
                        if strand_one == '+' and strand_two == '-':
                            contact_map[site_one - s, site_two - s, 0] += 1
                        """ OUTWARD """
                        if strand_one == '-' and strand_two == '+':
                            contact_map[site_one - s, site_two - s, 1] += 1
                        """ TANDEM P """
                        if strand_one == '+' and strand_two == '+':
                            contact_map[site_one - s, site_two - s, 2] += 1
                        """ TANDEM N """
                        if strand_one == '-' and strand_two == '-':
                            contact_map[site_one - s, site_two - s, 3] += 1
                # break will trigger both when site_one and site_two are both
                # ind mapped, and also when site_one is in one TSS (as defined
                # by chromatin_len) and site two is in another
                break

np.save('contact_map_' + str(slice_one) + '.npy', contact_map)