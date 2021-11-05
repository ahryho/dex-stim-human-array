
#!/bin/bash

dir="~/bio/code/mpip/dex-stim-human-array/data/snps/final_imputed_qc_snps"
dir="/binder/mgp/workspace/2020_DexStim_Array_Human/snps/imputed_qc/from_darina/ld_200_100"

# recode the Dex_genoData into .gen format + .sample
$dir/plink --bfile $dir/pca/ld_extracted/dex_geno_imputed_ld --recode oxford --out $dir/mtrx/dex_geno_imputed_ld
echo "recoded Dex_genoData"

# retrieve the SNP position - consisting of chromosome, ID and position
cut -f 2 -d' ' $dir/mtrx/dex_geno_imputed_ld.gen | sed 's/ /_/g' > $dir/mtrx/snps_positions.txt
echo "done with retrieving info about SNP postion"

# retrieve the info with three numbers describing one sample (1st: homozygote reference, 2nd: heterozygote, 3rd: homozygote alternative)
cut -f 6- -d' ' $dir/mtrx/dex_geno_imputed_ld.gen > $dir/mtrx/plinkgen_X.txt
echo "done with retrieving info about samples"

# Copy from slurngate to local 
# scp -C -r ahryhorzhevska@slurmgate:/binder/mgp/workspace/2020_DexStim_Array_Human/snps/imputed_qc/from_darina/ld_200_100/mtrx  ~/bio/code/mpip/dex-stim-human-array/data/snps/final_imputed_qc_snps