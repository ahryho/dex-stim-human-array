#!/bin/bash

if [ $# -ne 6 ]; then
    echo "$0: usage: Not enough arguments
          First argument: path to the R Script
          Second argument: result file name
	        Third argument: model
	        Fourth argument: partition
	        Fifth argiment: node
	        Sixth argument: memoty in Gb"
    exit 1
fi

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-differential-methyl-analysis
rslt_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/01_lme_models

beta_mtrx_fn=/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data_no_outlier/dex_methyl_beta_combat_mtrx.rds
pheno_fn=/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/pheno_full_for_kimono.csv

# out_fn=lmem_dex.out
# err_fn=.err

module load R

# job_name=lme_dex

r_script=$1
rslt_lmem_fn=$2
model=$3
partition=$4
node=$5
memory=$6

job_name=lmem_$3_dex
out_fn=lmem_dex_$3.out

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla $r_script $beta_mtrx_fn $pheno_fn $rslt_dir/$rslt_lmem_fn"
