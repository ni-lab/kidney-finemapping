#!/bin/bash

## Load modules and set environment variables:
module load python/3.8.8
module load cuda/10.0
source activate kidney_finemapping

export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
export PATH=$BASE_DIR/bin:$PATH
export PYTHONPATH=$BASE_DIR/bin:$PYTHONPATH

## Command(s) to run:
cd $BASE_DIR
merge_ism_summed_pos_shifts.py \
  /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts \
  -o /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts/all_chrs
