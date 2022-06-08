#!/bin/bash
#
#SBATCH --job-name=tca
#SBATCH --output=err_out/tca_%A_%a.out
#SBATCH --error=err_out/tca_%A_%a.err
#SBATCH --mem-per-cpu=100Gb 
#SBATCH --array=2-22%22    # Submit 20 tasks. Run max 20 concurrently
#SBATCH --part=pe,hp
#SBATCH --exclude=pe9

if [ $# -ne 2 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex)
          Second argument: path to the R Script"
    exit 1
fi

module load R

treatment=$1
r_script=$2

chr=$SLURM_ARRAY_TASK_ID

Rscript --vanilla $r_script $treatment $chr
