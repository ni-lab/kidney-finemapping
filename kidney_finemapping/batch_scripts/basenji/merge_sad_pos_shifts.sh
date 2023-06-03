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

#!/bin/bash
## Load modules and set environment variables:
export PATH=/clusterfs/nilah/richard/home/conda/envs/basenji_kidney_finemapping:$PATH
module load cuda/11.2

# source the conda.sh script:
source /global/software/sl-7.x86_64/modules/langs/python/3.7/etc/profile.d/conda.sh
conda activate basenji_kidney_finemapping
export BASE_DIR=/clusterfs/nilah/richard/home/kidney-finemapping

## Command(s) to run:
cd $BASE_DIR
kidney_finemapping/basenji/merge_sad.py out_dir/220513_variants/sad_shifts \
    -o out_dir/220513_variants/sad_shifts/all_chrs
