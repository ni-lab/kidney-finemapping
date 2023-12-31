#!/bin/bash
# Job name:
#SBATCH --job-name=preprocess_finemapped_variants
#
# Account:
#SBATCH --account=fc_nilah
#
# Partition:
#SBATCH --partition=savio2
#
# Number of nodes:
#SBATCH --nodes=1
#
# Number of tasks (one for each GPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task (please always specify the total number of processors twice the number of GPUs):
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=30:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
export PATH=/clusterfs/nilah/richard/home/conda/envs/basenji_kidney_finemapping:$PATH
module load cuda/11.2

# source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping

# Set base directory
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
cd $BASE_DIR

## Command(s) to run:
# 220620 variants
kidney_finemapping/basenji/preprocess_finemapped_variants.py \
    resources/data/220620_variants/raw/220620_tubule_centered_susie3_exonexclude_220617_FinemappedTubuleTest.txt \
    --variant_set 220620 \
    -o out_dir/220620_variants/susie/data/preprocessed

# 220513 variants
kidney_finemapping/basenji/preprocess_finemapped_variants.py \
    resources/data/220513_variants/raw/220513_gwas_replicating_loci_top_variants_susie5_finemap5_01_hg38.txt \
    --variant_set 220513 \
    -o out_dir/220513_variants/data/preprocessed
