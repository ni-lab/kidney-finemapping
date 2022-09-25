#!/bin/bash
# Job name:
#SBATCH --job-name=basenji_sad_pos_shifts
#
# Account:
#SBATCH --account=co_nilah
#
# Partition:
#SBATCH --partition=savio2_gpu
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
#SBATCH --time=12:00:00
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
source activate basenji
export BASENJIDIR=/global/home/users/richardshuai/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH

## Command(s) to run:
cd /clusterfs/nilah/richard/home/basenji
for CHROM in {1..22}
do
basenji_sad_pos_shifts.py --stats SAD,REF,ALT -f /clusterfs/nilah/richard/genomes/hg38.ml.fa -o /clusterfs/nilah/richard/kidney_data/220513_variants/sad_pos_shifts/chr${CHROM} --rc --shifts "1,0,-1" -t /clusterfs/nilah/richard/kidney_data/targets/kidney_sc_wigs_hg38.txt /clusterfs/nilah/richard/kidney_data/models/params_sc_kidney_regression.json /clusterfs/nilah/pooja/kidney_data/train/regression_chr_models/train_bigwigs_chr${CHROM}/model_best.h5  /clusterfs/nilah/richard/kidney_data/220513_variants/data/processed/snps_by_chrom/chr${CHROM}_snps.vcf
done
