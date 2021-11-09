#!/bin/bash

if [ $# -ne 6 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex, delta)
          Second argument: source folder path 
	        Third argument: result folder path
	        Fourth argument: partition
	        Fifth argument: node
	        Sixth argument: memory in Gb"
    exit 1
fi

module load R

treatment=$1
src_dir=$2
rslt_dir=$3
partition=$4
node=$5
memory=$6

job_name=eqtl_$1
out_fn=eqtl_$1.out

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/code/integrative/eqtl/03_run_eqtl.R $src_dir $rslt_dir"

# ./03_run_eqtl.sh dex /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/ /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/ pe 7 500Gb