#!/bin/bash
## Load modules and set environment variables:
module load python
module load cuda/10.0
source activate basenji
export BASENJIDIR=/clusterfs/nilah/richard/home/basenji/
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH

## Command(s) to run:
cd /clusterfs/nilah/richard/home/basenji
merge_ism_summed_pos_shifts.py /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts -o /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts/all_chrs 
