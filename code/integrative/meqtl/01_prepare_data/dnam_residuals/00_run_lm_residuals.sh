#!/bin/bash

if [ $# -ne 6 ]; then
    echo "$0: usage: Not enough arguments
          First argument: path to the R Script
          Second argument: result file name
	        Third argument: treatment (dex or veh)
	        Fourth argument: partition
	        Fifth argiment: node
	        Sixth argument: memoty in Gb"
    exit 1
fi

beta_mtrx_fn="/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_"
pheno_fn="/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"

module load R

# job_name=lme_dex

r_script=$1
rslt_lmem_fn=$2
treatment=$3
partition=$4
node=$5
memory=$6

job_name=resid_$3
out_fn=resid_$3.out

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla $r_script $beta_mtrx_fn $pheno_fn $rslt_dir/$rslt_lmem_fn $treatment"

# ./00_run_lm_residuals.sh lm_residuals_svs_pcs_smoke.R /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals_dex dex pe 7 900Gb

