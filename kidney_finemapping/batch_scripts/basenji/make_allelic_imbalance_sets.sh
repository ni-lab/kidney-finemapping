#!/bin/bash
# Job name:
#SBATCH --job-name=make_allelic_imbalance_sets
#
# Account:
#SBATCH --account=ac_nilah
#
# Partition:
#SBATCH --partition=savio3
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
module load cuda/11.2

# source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping

# Set base directory
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
cd $BASE_DIR

## Command(s) to run:

# Make allelic imbalance sets for LOH, PT, and DT
NEG_MULT=7
THRESH=0.01
for TARGET in LOH PT DT; do
    echo "Processing ${TARGET}..."
    kidney_finemapping/basenji/make_allelic_imbalance_sets.py \
        out_dir/allelic_imbalance/data/raw/astestq10tab/all_${TARGET}q10.tsv \
        out_dir/sc_atac_seq/${TARGET}_peaks.narrowPeak \
        --neg_mult ${NEG_MULT} \
        --n_bins 20 \
        --thresh ${THRESH} \
        -o out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}
done

# Make allelic imbalance sets for combined, Tubule
NEG_MULT=2
for TARGET in combined Tubule; do
    echo "Processing ${TARGET}..."
    kidney_finemapping/basenji/make_allelic_imbalance_sets.py \
        out_dir/allelic_imbalance/data/raw/astestq10tab/all_${TARGET}q10.tsv \
        out_dir/sc_atac_seq \
        --tubule_peaks \
        --neg_mult ${NEG_MULT} \
        --n_bins 20 \
        --thresh ${THRESH} \
         -o out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}
done
