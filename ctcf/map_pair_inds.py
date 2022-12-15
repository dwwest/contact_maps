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
with open("intersected_GSM2418860_expanded-pubWTmESC-mapped-filt.pairs") as f:
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
with open("intersected_GSM2418860_expanded.bed") as f:
    for ind,line in enumerate(f):
        if ind > 0:
            current = line.strip('\n').split('\t')
            if current[0] not in ctcf.keys():
                ctcf[current[0]] = [int(current[1])]
            else:
                ctcf[current[0]].append(int(current[1]))
print('Sites loaded...')

current_chrom = ''
wrap = 120
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
        # find which ctcf site from granges bed that the site aligns to
        for s in ctcf[chrom_one]:
            # check to make sure we are only looking at mononucleosomes
            if abs(site_two - site_one) < 400:
                # check if the first and second sites are both within the expanded ctcf area
                if site_one >= s and site_one <= chromatin_len + s:
                    if site_two >= s and site_two <= chromatin_len + s:
                        # for in, out, and tandem, 
                        """ INWARD """
                        if strand_one == '+' and strand_two == '-':
                            if chromatin_len >= site_one - s + wrap >= 0 and chromatin_len >= site_two - s - wrap >= 0:
                                contact_map[site_one - s + wrap, site_two - s - wrap, 0] += 1
                        """ OUTWARD """
                        if strand_one == '-' and strand_two == '+':
                            if chromatin_len >= site_one - s - wrap >= 0 and chromatin_len >= site_two - s + wrap >= 0:
                                contact_map[site_one - s - wrap, site_two - s + wrap, 1] += 1
                        """ TANDEM P """
                        if strand_one == '+' and strand_two == '+':
                            if chromatin_len >= site_one - s + wrap >= 0 and chromatin_len >= site_two - s + wrap >= 0:
                                contact_map[site_one - s + wrap, site_two - s + wrap, 2] += 1
                        """ TANDEM N """
                        if strand_one == '-' and strand_two == '-':
                            if chromatin_len >= site_one - s - wrap >= 0 and chromatin_len >= site_two - s - wrap >= 0:
                                contact_map[site_one - s - wrap, site_two - s - wrap, 3] += 1
                    # break will trigger both when site_one and site_two are both
                    # ind mapped, and also when site_one is in one TSS (as defined
                    # by chromatin_len) and site two is in another
                    break

to_flip = np.zeros([chromatin_len + 1, chromatin_len + 1])
for i in range(chromatin_len):
    for j in range(i+1, chromatin_len+1):
        to_flip[i, j] = 1

to_save = np.zeros(np.shape(contact_map))
for i in range(4):
    to_add = np.multiply(to_flip, contact_map[:,:,i])
    to_save[:,:,i] = np.flip(np.flip(to_add, 1), 0) + contact_map[:,:,i]
    flattening = np.flip(np.flip(to_flip, 1), 0)
    to_save[:,:,i] = np.multiply(flattening, to_save[:,:,i])
np.save('contact_map_' + str(slice_one) + '.npy', to_save)