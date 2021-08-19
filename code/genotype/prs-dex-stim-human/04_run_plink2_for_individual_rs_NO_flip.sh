#!/bin/bash
#
#SBATCH --job-name=prs-ind-rs
#SBATCH --output=err_out/prs_individual_%A_%a.out
#SBATCH --error=err_out/prs_individual_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mem=5000
#SBATCH --exclude=hp01,hp02

# Complete path to the target data (bim/bed files)
cohort="/binder/common/genotypes/qc_imputations/DexStim_Mar2015/bggt/dex_imputed_bggt_all_qc"
# List of cohorts
cohort_short="MARS_DexStim"
ncohort=${#cohort[@]}

# Path to the folder with PRS-CS output
prs_res_dir="/binder/mgp/datasets/2020_DexStim_Array_Human/PRS"

# List of diseases/traits that are in the folder containing PRS-CS results
# traits=( $(ls -F $prs_res_dir | grep / | sed 's#/##') )
traits=("AD" "ADHD")
ntraits=${#traits[@]}

# Get trait index for each job id
itrait=$(($SLURM_ARRAY_TASK_ID%$ntraits))

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "PRS result dir: " $prs_res_dir
echo "Cohort: " $cohort_short
echo "Cohort path to bim/bed files: " $cohort
echo "Trait: " "${traits[$itrait]}"
echo "Trait index: " $itrait

# Run PLINK2 to calculate overall risk of each individual in cohorts for the different diseases/traits
/binder/mgp/code/plink2/plink2 --score "${prs_res_dir}/${traits[$itrait]}/${cohort_short}_${traits[$itrait]}_all_chromosomes.txt" 2 4 6 \
                               --out "${prs_res_dir}/${traits[$itrait]}/${cohort_short}_${traits[$itrait]}_individual_rs_plink_out_no_flip" \
                               --bfile $cohort

# Remove concatenated file
# rm "${prs_res_dir}/${traits[$itrait]}/${cohort_short}_${traits[$itrait]}_all_chromosomes.txt"
