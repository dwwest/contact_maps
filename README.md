# contact_maps

My code to make cutsom ctcf and tss contact maps from micro-c data

For the CTCF Pipeline:
1) First run text_to_bed.py to make bed file from peaks file for CTCF occupancy
2) Next, run expand_granges.py to expand the granges of ctcf out by 800 bp to 
    either side of the CTCF site -- we want to capture structure of chromatin fiber
    on both sides, not just the CTCF site itself.
3) Run occupied-intersect.sh to intersect the expanded granges with the bed
    peaks to get a set of occupied granges.
4) Run anno-intersect.sh on the expanded, occupied CTCF sites to intersect with 
    the Micro-C pairs and bam files to get Micro-C signals at occupied CTCF sites.
5) Run map_pair_inds.sh which runs several parallel instances of map_pair_inds.py 
    to align all pairs with both ends inside the ctcf site and determine
    sum of pairs for each bp and the orientations.  Account for wrapping in order
    to get, not ends of reads, but positions of entry and exit base pairs as given
    by whether reads are in, out, or tandem.
6) Run make_heatmap in order to generate bp resolution heatmap, sum bp resolution to get
    nucleosome resolution heatmap, and pileup histograms.