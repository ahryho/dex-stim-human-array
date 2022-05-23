#!/bin/bash
#
#SBATCH --job-name=tca
#SBATCH --output=err_out/tca_%A_%a.out
#SBATCH --error=err_out/tca_%A_%a.err
#SBATCH --array=1-20%20    # Submit 20 tasks. Run max 20 concurrently
#SBATCH --part=pe
#SBATCH --exclude=pe6

if [ $# -ne 3 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex)
          Second argument: path to the R Script"
    exit 1
fi

module load R

treatment=$1
r_script=$2

chr=$SLURM_ARRAY_TASK_ID

sbatch --wrap="Rscript --vanilla $r_script $treatment $chr"
