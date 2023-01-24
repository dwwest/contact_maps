import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import glob
import re
import sys
import os
from scipy.signal import butter
from scipy.signal import lfilter

to_save_dir = '../../wlcsim/isd_input/micro-c_data/hsieh_ctcf/'
cm_files = glob.glob('contact_map_*.npy')
chrom_length = 1600

cm = np.zeros([chrom_length + 1, chrom_length + 1, 4])
for f in cm_files:
    current_cm = np.load(f)
    cm = cm + current_cm

to_plot = np.zeros([chrom_length + 1, chrom_length + 1, 3])
to_plot[:,:,0:2] = cm[:,:,0:2]
to_plot[:,:,2] = cm[:,:,2] + cm[:,:,3]
max_contacts = np.max(to_plot)
to_plot = 1-(to_plot/max_contacts) # get the secondary colors
for i in range(3):
    to_plot[:,:,i] = np.transpose(to_plot[:,:,i])
# CYAN = IN
# YELLOW = OUT
# MAGENTA = TANDEM

plt.imshow(to_plot)
plt.savefig('heatmap.pdf')
plt.show()

to_bin = cm[240:1361,240:1361,:]
bin_boundaries = [0, 160, 320, 530, 590, 800, 960, 1120]
ctcf_ind = 3
binned_array = np.zeros([36, 4])
ind = 0
for inda in range(len(bin_boundaries)-1):
    for indb in range(len(bin_boundaries)-1):
        # leave out the bin where the ctcf site is
        if not inda == 3 and not indb == 3:
            for orient in range(4):
                binned_array[ind, orient] = np.sum(to_bin[bin_boundaries[inda]:bin_boundaries[inda+1],\
                    bin_boundaries[indb]:bin_boundaries[indb+1],orient])
            ind += 1
binned_array = np.reshape(binned_array, [6, 6, 4])
for i in range(4):
    binned_array[:,:,i] = binned_array[:,:,i].T
to_plot_binned = np.zeros([6, 6, 3])
to_plot_binned[:,:,0:2] = binned_array[:,:,0:2]
to_plot_binned[:,:,2] = binned_array[:,:,2] + binned_array[:,:,3]
max_contacts = np.max(to_plot_binned)
to_plot_binned = 1-(to_plot_binned/max_contacts)
plt.imshow(to_plot_binned)
plt.savefig('binned_heatmap.pdf')
plt.show()

if not os.path.isdir(to_save_dir):
    os.mkdir(to_save_dir)
with open(to_save_dir + '/data.txt', 'w+') as f:
    for k in range(4):
        for i in range(len(binned_array)-1):
            for j in range(i+1, len(binned_array)):
                print(i, j)
                f.write(str(binned_array[i, j, k]) + '\n')

freqs_by_pos = np.zeros(chrom_length+1)
for i, row in enumerate(cm):
    for j, pos in enumerate(row):
        freqs_by_pos[i] += np.sum(pos)
        freqs_by_pos[j] += np.sum(pos)

# lowpass filter
cutoff = 0.1
nyq = 0.5 * 1 # for 1bp
normal_cutoff = cutoff / nyq
b, a = butter(2, normal_cutoff, btype='low', analog=False)
y = lfilter(b, a, freqs_by_pos)

xs = [i for i in range(-801, 800)]
plt.plot(xs, freqs_by_pos)
plt.plot(xs, y)
plt.vlines([-577, -544, -417, -384, -257, -227, -100, 100, 227, 257, 384, 417, 544, 577], [0, 0, 0], [80000, 80000, 80000], colors = ['r', 'r', 'r'])
plt.xlabel('Fragment end position')
plt.ylabel('Frequency')
plt.title('CTCF Micro-C nucleosome positions')
plt.savefig('ctcf_freqs_peaks.pdf')
plt.show()

xs = [i for i in range(-801, 800)]
plt.plot(xs, freqs_by_pos)
plt.plot(xs, y)
plt.vlines([-560, -400, -240, -30, 30, 240, 400, 560], [0, 0, 0], [80000, 80000, 80000], colors = ['r', 'r', 'r'])
plt.xlabel('Fragment end position')
plt.ylabel('Frequency')
plt.title('CTCF Micro-C fragment end bin bounds')
plt.savefig('ctcf_freqs_bins.pdf')
plt.show()
