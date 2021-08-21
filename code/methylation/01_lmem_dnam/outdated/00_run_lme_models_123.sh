#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-differential-methyl-analysis
rslt_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models

input_param=$src_dir/input_parameters.csv

rslt_lme=$rslt_dir/lme_dex

# out_fn=.out
# err_fn=.err

# module load R

sbatch --job-name=lme_dex_sex --part=pe --mem=300G --output=${rslt_lme}_sex.out --error=${rslt_lme}_sex.err --wrap="Rscript --vanilla $src_dir/03_lme_models/01_lme_model_sex.R $input_param ${rslt_lme}_sex.txt" &
sbatch --job-name=lme_dex_all --part=pe --mem=300G --output=${rslt_lme}_all.out --error=${rslt_lme}_all.err --wrap="Rscript --vanilla $src_dir/03_lme_models/02_lme_model_all.R $input_param ${rslt_lme}_all.txt" &
sbatch --job-name=lme_dex_all_cells --part=pe --mem=300G --output=${rslt_lme}_all_plus_cell_counts.out --error=${rslt_lme}_all_plus_cell_counts.err --wrap="Rscript --vanilla $src_dir/03_lme_models/03_lme_model_all_plus_cell_counts.R $input_param ${rslt_lme}_all_plus_cell_counts.txt" &

wait
