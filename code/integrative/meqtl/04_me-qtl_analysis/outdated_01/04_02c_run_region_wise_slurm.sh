#!/bin/bash

treatment="veh"

for i in 0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 110000 120000
do
  sbatch ./04_02c_get_independent_me-qtl_snps_region_wise_slurm.sh $i $treatment
done

find ${treatment}/. -type f -name "*.txt" -exec cat {} + > me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

