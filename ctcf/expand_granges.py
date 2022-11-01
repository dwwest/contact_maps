import numpy as np

"""
This is here to take the granges file for ctcf sites and extend the grange out by 800 bp on either side
so that when I run anno-intersect.sh, I can intersect for micro-c pairs that fall not directly on the 
ctcf site but on +2 and +3 nucleosomes as well
"""

bp_to_extend = 790 # how many bp to extend out from granges
with open('granges.bed', 'r') as f_og, open('granges_expanded.bed', 'w+') as f_expand:
    f_expand.write('# \n')
    for ind,line in enumerate(f_og):
        if ind > 0:
            current = line.strip('\n').split('\t')
            # extend range back by bp_to_extend
            if int(current[1]) - bp_to_extend >= 0:
                current[1] = str(int(current[1]) - bp_to_extend)
            # extend range forward by bp_to_extend
            current[2] = str(int(current[2]) + bp_to_extend)
            # write larger grange to granges_expanded file
            f_expand.write('\t'.join(current[0:3]) + '\n')
