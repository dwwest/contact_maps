import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import glob
import re
import sys
import os

input_dir = '/rugpfs/fs0/risc_lab/scratch/dwest/isd_wlcsim/workspace/isd_rise_per_nuc_micro-c/wlcsim/input/'
lls = [33., 33., 33., 200., 33., 33., 33.]
discs = [i/4 for i in lls]

with open(input_dir + 'ctcf', 'w+') as f:
    for i in range(len(lls)):
        f.write(str(lls[i]) + '\t' + str(discs[i]) + '\n')

