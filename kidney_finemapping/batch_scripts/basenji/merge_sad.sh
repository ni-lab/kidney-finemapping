#!/bin/bash
## Load modules and set environment variables:
export PATH=/clusterfs/nilah/richard/home/conda/envs/basenji_kidney_finemapping:$PATH
module load cuda/11.2

# Source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping

# Set base directory
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
cd $BASE_DIR

## Command(s) to run:
kidney_finemapping/basenji/merge_sad.py out_dir/220620_variants/susie/sad/ \
  -n 22 \
  --vcf \
  -o out_dir/220620_variants/susie/sad/all_chrs
