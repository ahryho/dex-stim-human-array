#!/bin/bash

#SBATCH --job-name=prs-flip-prs-output
#SBATCH --output=err_out/prs_flip_%A_%a.out
#SBATCH --error=err_out/prs_flip_%A_%a.err
#SBATCH --array=0-21
#SBATCH --mem=5000
#SBATCH --exclude=hp01,hp02

# Path to the folder with PRS-CS output
prs_res_dir="/binder/mgp/datasets/2020_DexStim_Array_Human/PRS"

# List of cohorts
cohort_short=("MARS_DexStim")
ncohort=${#cohort_short[@]}

# List of diseases/traits that are in the folder containing PRS-CS results
# traits=( $(ls -F $prs_res_dir | grep / | sed 's#/##') )
traits=("AD" "ADHD")
ntraits=${#traits[@]}

# Get cohort and trait index for each job id
icohort=$((SLURM_ARRAY_TASK_ID / ntraits))
icohort=$icohort|cut -f1 -d"."
itrait=$(($SLURM_ARRAY_TASK_ID%$ntraits))

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "PRS result dir: " $prs_res_dir
echo "Cohort: " "${cohort_short[$icohort]}"
echo "Trait: " "${traits[$itrait]}"
echo "Cohort index: " $icohort
echo "Trait index: " $itrait

# Flip effect sizes in PRS-CS output if negative
Rscript --vanilla 03_flip_prs_output.R "${traits[$itrait]}" "${cohort_short[$icohort]}" $prs_res_dir
