#!/bin/bash
#
#SBATCH --job-name=PRS
#SBATCH --output=/home/ahryhorzhevska/mpip/code/prs-dex-stim-human/err_out/prs_%A_%a.out
#SBATCH --error=/home/ahryhorzhevska/mpip/code/prs-dex-stim-human/err_out/prs_%A_%a.err
#SBATCH --array=0
#SBATCH --mem=12000
#SBATCH --exclude=hp01,hp02

# Complete path to the target data (bim file)
cohort="/binder/common/genotypes/qc_imputations/DexStim_Mar2015/bggt/dex_imputed_bggt_all_qc"
# cohort_short="MARS_DexStim_Mar2015"
ncohort=${#cohort[@]}

# Complete path to the GWAS statistics
gwas_summaries="/binder/mgp/datasets/2020_PRS_Covid19/GWAS_sumstats/GWAS_summaries.tsv"

# List of diseases/traits and the respective GWAS sample size taken from summary file, remove header line
traits=("ASD")
ntraits=( ${#traits[@]} )
IFS=$'\n'
nsamples=( $(cut -f 3,8 $gwas_summaries |  grep -- "${traits[*]}" | cut -f1) )
phicode=( $(cut -f 8,10 $gwas_summaries |  grep -- "${traits[*]}" | cut -f2) )

# Get a trait index for each job id
itrait=$(($SLURM_ARRAY_TASK_ID%$ntraits))

# Create directory for storing results of each cohort
prs_res_dir="/binder/mgp/datasets/2020_DexStim_Array_Human/PRS/${traits[$itrait]}"
mkdir -p ${prs_res_dir}

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "Cohort: " $cohort
echo "Trait: " "${traits[$itrait]}"
echo "Trait index: "  $itrait
echo "Number of samples: " "${nsamples[$itrait]}"
echo "Phicode: " "${phicode[$itrait]}"
echo "Result directory: " ${prs_res_dir}

# Run PRS-CS
gwas_sumstats_prefix="/binder/mgp/datasets/2020_PRS_Covid19/GWAS_sumstats"
if [ ${phicode[$itrait]} == 2 ]
then
        python3 /binder/common/PRS-CS/PRScs.py \
            --ref_dir=/binder/common/PRS-CS/ldblk_1kg_eur/ \
            --bim_prefix=$cohort \
            --sst_file="${gwas_sumstats_prefix}/${traits[$itrait]}/${traits[$itrait]}.flipRisk" \
            --n_gwas=${nsamples[$itrait]} \
            --out_dir="${prs_res_dir}/${traits[$itrait]}"
elif [ ${phicode[$itrait]} == 1 ]
then
        python3 /binder/common/PRS-CS/PRScs.py \
            --ref_dir=/binder/common/PRS-CS/ldblk_1kg_eur/ \
            --bim_prefix=$cohort \
            --sst_file="${gwas_sumstats_prefix}/${traits[$itrait]}/${traits[$itrait]}.flipRisk" \
            --n_gwas=${nsamples[$itrait]} \
            --phi=1e-2 \
            --out_dir="${prs_res_dir}/${traits[$itrait]}"
else
        python3 /binder/common/PRS-CS/PRScs.py \
            --ref_dir=/binder/common/PRS-CS/ldblk_1kg_eur/ \
            --bim_prefix=$cohort \
            --sst_file="${gwas_sumstats_prefix}/${traits[$itrait]}/${traits[$itrait]}.flipRisk" \
            --n_gwas=${nsamples[$itrait]} \
            --phi=1e-4 \
            --out_dir="${prs_res_dir}/${traits[$itrait]}"
fi
