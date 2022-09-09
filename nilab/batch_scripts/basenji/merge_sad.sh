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
merge_sad.py -n 22 -o /clusterfs/nilah/richard/kidney_data/220620_variants/susie/sad/all_chrs --stats REF,ALT,SAD --vcf /clusterfs/nilah/richard/kidney_data/220620_variants/susie/sad/
