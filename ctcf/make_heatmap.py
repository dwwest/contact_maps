import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import glob
import re
import sys

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

