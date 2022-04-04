#!/bin/bash

if [ $# -ne 7 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex, delta)
          Second argument: source folder path 
	        Third argument: result folder path
	        Fourth argument: partition
	        Fifth argument: node
	        Sixth argument: memory in Gb
	        Seventh argument: path to the R Script"
    exit 1
fi

module load R

treatment=$1
src_dir=$2
rslt_dir=$3
partition=$4
node=$5
memory=$6
r_script=$7

job_name=eqtm_$1
out_fn=$job_name.out

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla $r_script $treatment $src_dir $rslt_dir"

# ./03_run_eqtm.sh dex /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/ /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/ pe 7 400Gb 03_run_eqtm.R
