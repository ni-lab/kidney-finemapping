#!/bin/bash
## Load modules and set environment variables:

## Command(s) to run:
module load python/3.8.8
module load cuda/10.0
conda activate kidney_finemapping

export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
export PATH=$BASE_DIR/bin:$PATH
export PYTHONPATH=$BASE_DIR/bin:$PYTHONPATH

## Command(s) to run:
cd $BASE_DIR
merge_sad.py /clusterfs/nilah/richard/kidney_data/220620_variants/susie/sad/ \
  -n 22 \
  -o /clusterfs/nilah/richard/kidney_data/220620_variants/susie/sad/all_chrs
