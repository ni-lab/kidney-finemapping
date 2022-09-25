#!/bin/bash
# Job name:
#SBATCH --job-name=basenji_merge_sad_pos_shifts
#
# Account:
#SBATCH --account=fc_nilah
#
# Partition:
#SBATCH --partition=savio
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
#SBATCH --time=1:00:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
module load python/3.8.8
module load cuda/10.0
source activate basenji
export BASENJIDIR=/global/home/users/richardshuai/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH

## Command(s) to run:
merge_sad.py -o /clusterfs/nilah/richard/kidney_data/220513_variants/sad_pos_shifts/all_chrs --stats REF,ALT,SAD /clusterfs/nilah/richard/kidney_data/220513_variants/sad_pos_shifts/
