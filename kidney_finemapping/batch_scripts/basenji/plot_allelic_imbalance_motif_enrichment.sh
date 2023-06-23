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
#SBATCH --time=24:00:00
#
# Send stdout and stderr to logs
#SBATCH --output=log_files/job_%j.out
#SBATCH --error=log_files/job_%j.err
#

## Write batch script to logs
scontrol write batch_script $SLURM_JOB_ID log_files/job_$SLURM_JOB_ID.sh

## Load modules and set environment variables:
export PATH=/clusterfs/nilah/richard/home/conda/envs/basenji_kidney_finemapping:$PATH

# source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping

# Set base directory
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping
cd $BASE_DIR

for TARGET in PT LOH DT; do
  python3 kidney_finemapping/basenji/plot/plot_allelic_imbalance_motif_enrichment.py \
  resources/data/tf_pseudobulk_Pseudobulk_Wilson_TF_analysis.csv \
  out_dir/allelic_imbalance/motif_enrichment/${TARGET}_variants_neg7x_q0.01/hypergeom_per_motif.tsv \
  --cell_type ${TARGET} \
  -o out_dir/allelic_imbalance/plots/motif_enrichment/${TARGET}_variants
done
