#!/bin/sh

#SBATCH --job-name=prs-concat
#SBATCH --output=err_out/prs_concat_%A_%a.out
#SBATCH --error=err_out/prs_concat_%A_%a.err
#SBATCH --array=0
#SBATCH --mem=5000
#SBATCH --exclude=hp01,hp02

prs_res_dir="/binder/mgp/datasets/2020_DexStim_Array_Human/PRS"
cohort_short="MARS_DexStim"

# List of diseases/traits that are in the folder containing PRS-CS results
# traits=( $(ls -F $prs_res_dir | grep / | sed 's#/##') )
traits=("AD")
ntraits=${#traits[@]}

# Get trait index for each job id
itrait=$(($SLURM_ARRAY_TASK_ID%$ntraits))

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "PRS result dir: " $prs_res_dir
echo "Cohort: " $cohort_short
echo "Trait: " "${traits[$itrait]}"
echo "Trait index: " $itrait

# Concatenate chromosome summary statistics files to one file
cat ${prs_res_dir}/${traits[$itrait]}/${traits[$itrait]}*.txt > ${prs_res_dir}/${traits[$itrait]}/${cohort_short}_${traits[$itrait]}_all_chromosomes.txt

# Move single chrom into folder with cohort name
mkdir -p ${prs_res_dir}/${traits[$itrait]}/${cohort_short}
mv ${prs_res_dir}/${traits[$itrait]}/${traits[$itrait]}*.txt ${prs_res_dir}/${traits[$itrait]}/${cohort_short}
