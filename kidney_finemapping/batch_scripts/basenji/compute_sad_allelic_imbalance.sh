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
NEG_MULT=7
THRESH=0.01
for TARGET in LOH PT DT; do
  for SET in pos neg; do
    for CHR in {1..22}; do
      CHROM=chr${CHR}
      kidney_finemapping/basenji/compute_sad.py \
        resources/model_params/params_sc_kidney_regression.json \
        resources/models/train_bigwigs_${CHROM}/model_best.h5 \
        out_dir/allelic_imbalance/data/preprocessed/${TARGET}_variants_neg${NEG_MULT}x_q${THRESH}/snps_by_chrom/${SET}_set_${CHROM}.vcf \
        -f resources/genomes/hg38.ml.fa \
        --rc \
        --shifts "1,0,-1" \
        -t resources/targets/kidney_sc_wigs_hg38.txt \
        -o out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${Q_THRESH}/${SET}_sad/${CHROM}
    done
    # Merge SAD files across chromosomes
    kidney_finemapping/basenji/merge_sad.py \
      out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${Q_THRESH}/${SET}_sad \
      --vcf \
      -o out_dir/allelic_imbalance/sad/${TARGET}_neg${NEG_MULT}x_q${Q_THRESH}/${SET}_sad/all_chrs
  done
done

# Merge SAD files for each chr
for TARGET in "${TARGETS[@]}"; do
        for SET in "${SETS[@]}"; do

        done
done

