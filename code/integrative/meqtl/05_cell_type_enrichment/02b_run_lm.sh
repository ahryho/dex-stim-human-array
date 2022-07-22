#!/bin/bash

if [ $# -ne 4 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex)
          Second argument: partition
	        Third argument: node
	        Forth argument: memory in Gb"
    exit 1
fi

module load R

treatment=$1
partition=$2
node=$3
memory=$4

job_name=cell_type_lm_$1
out_fn=cell_type_lm_$1.out

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla /home/ahryhorzhevska/mpip/projects/dex-stim-human-array/code/integrative/meqtl/05_cell_type_enrichment/02b_lm.R $treatment"