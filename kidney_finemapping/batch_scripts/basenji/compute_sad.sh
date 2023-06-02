#!/bin/bash
# Job name:
#SBATCH --job-name=compute_sad
#
# Account:
#SBATCH --account=ac_nilah
#
# Partition:
#SBATCH --partition=savio2_1080ti
#
# Number of nodes:
#SBATCH --nodes=1
#
# Number of tasks (one for each GPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task (please always specify the total number of processors twice the number of GPUs):
#SBATCH --cpus-per-task=2
#
#Number of GPUs, this can be in the format of "gpu:[1-4]", or "gpu:K80:[1-4] with the type included
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
export PATH=/clusterfs/nilah/richard/home/conda/envs/basenji_kidney_finemapping:$PATH
export TF_CPP_MIN_LOG_LEVEL=2
module load cuda/11.2

# Source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping

# Set base directory
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
cd $BASE_DIR

## Command(s) to run:
CHROM=chrXY
kidney_finemapping/basenji/compute_sad.py \
  /clusterfs/nilah/richard/kidney_data/models/params_sc_kidney_regression.json \
  /clusterfs/nilah/pooja/kidney_data/train/regression_chr_models/train_bigwigs_${CHROM}/model_best.h5 \
  out_dir/220620_variants/susie/data/preprocessed/snps_by_chrom/${CHROM}_snps.vcf \
  -f /clusterfs/nilah/richard/genomes/hg38.ml.fa \
  --rc \
  --shifts "1,0,-1" \
  -t /clusterfs/nilah/richard/kidney_data/targets/kidney_sc_wigs_hg38.txt \
  -o out_dir/220620_variants/susie/sad/${CHROM}
