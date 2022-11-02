#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1

source activate wlcsim

task_id=$SLURM_ARRAY_TASK_ID
task_count=$SLURM_ARRAY_TASK_COUNT

num_lines=`wc -l < intersected_GSM2418860_expanded-pubWTmESC-mapped-filt.pairs`
slice=$((num_lines/task_count*task_id))
last_id=$((task_id-1))
prev_slice=$((num_lines/task_count*last_id))
last=false
if [ "$task_id" == "$task_count" ]; then
    last=true
fi
python map_pair_inds.py $prev_slice $slice $last