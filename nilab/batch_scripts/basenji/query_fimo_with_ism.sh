#!/bin/bash
# Job name:
#SBATCH --job-name=basenji_fimo_with_ism
#
# Account:
#SBATCH --account=ac_nilah
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
#SBATCH --time=2:00:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
module load python
module load cuda/10.0
source activate basenji
module load meme
export BASENJIDIR=/clusterfs/nilah/richard/home/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH

## Command(s) to run:
cd /clusterfs/nilah/richard/home/basenji

MIN_THRESH=1e-3
RATIO_THRESH=0.5
query_fimo_with_ism.py --min_thresh ${MIN_THRESH} --ratio_thresh ${RATIO_THRESH} --stats max -t /clusterfs/nilah/richard/kidney_data/targets/kidney_sc_wigs_hg38.txt -o /clusterfs/nilah/richard/kidney_data/220513_variants/fimo/fimo_scored_by_ism_combined_dbs_${MIN_THRESH}_${RATIO_THRESH} /clusterfs/nilah/richard/motif_databases/databases_for_ism/meme_combined/meme_combined.meme /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts/all_chrs /clusterfs/nilah/richard/kidney_data/220513_variants/sad/all_chrs/sad.h5