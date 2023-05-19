#!/bin/bash
# Job name:
#SBATCH --job-name=basenji_ism_summed_pos_shifts
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
#SBATCH --gres=gpu:1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
conda init bash
source ~/.bashrc
module load python/3.8.8
module load cuda/10.0
conda activate kidney_finemapping

export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
export PATH=$BASE_DIR/bin:$PATH
export PYTHONPATH=$BASE_DIR/bin:$PYTHONPATH

## Command(s) to run:
CHROM=chr2
cd $BASE_DIR
kidney_finemapping/basenji/ism_summed_pos_shifts.py \
  /clusterfs/nilah/richard/kidney_data/models/params_sc_kidney_regression.json \
  /clusterfs/nilah/pooja/kidney_data/train/regression_chr_models/train_bigwigs_${CHROM}/model_best.h5 \
  /clusterfs/nilah/richard/kidney_data/220513_variants/data/processed/snps_by_chrom/${CHROM}_snps.vcf \
  -f /clusterfs/nilah/richard/genomes/hg38.ml.fa \
  --rc \
  --shifts "1,0,-1" \
  -t /clusterfs/nilah/richard/kidney_data/targets/kidney_sc_wigs_hg38.txt \
  -o /clusterfs/nilah/richard/refactor/kidney_data/220513_variants/ism_summed_pos_shifts/${CHROM}
