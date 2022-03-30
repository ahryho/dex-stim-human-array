#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono
data_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/eQTM/

# module load R

partition=hp
node=02
memory=100G
job_name=eqtm_

# treatment="veh"
# sbatch --job-name=$job_name$treatment --part=$partition --nodelist=$partition$node --mem=$memory \
# 	--wrap="Rscript --vanilla $src_dir/04_eqtm/03_run_eqtm.R $treatment $data_dir"

treatment="dex"
sbatch --job-name=$job_name$treatment --part=$partition --nodelist=$partition$node --mem=$memory \
	--wrap="Rscript --vanilla $src_dir/04_eqtm/03_run_eqtm.R $treatment $data_dir"
