#!/bin/bash
# Job name:
#SBATCH --job-name=basenji_ism_summed_pos_shifts
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
#SBATCH --cpus-per-task=2
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
source activate basenji

export BASENJIDIR=/clusterfs/nilah/richard/home/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH

## Command(s) to run:
cd /clusterfs/nilah/richard/home/basenji
basenji_plot_ism.py /clusterfs/nilah/richard/kidney_data/220513_variants/ism_summed_pos_shifts/all_chrs /clusterfs/nilah/richard/kidney_data/220513_variants/sad/all_chrs/sad.h5 -l 20 -g -t resources/targets/kidney_sc_wigs_hg38.txt -o /clusterfs/nilah/richard/kidney_data/220513_variants/ism_plots
