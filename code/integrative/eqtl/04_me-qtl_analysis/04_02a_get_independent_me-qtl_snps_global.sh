cd /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL

snp_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/snps

treatment=$1

echo "SNP P" > ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

awk 'NR != 1 {print $1 "\t" $5}' me-qtl_cis_result_${treatment}_fdr_005.csv | sort | uniq >> ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

# sed -i '1 i\SNP\tP' ind_me-qtl_cis/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

plink --bfile $snp_dir/Dex_genoData_SNPs --clump ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv --clump-r2 0.2 --clump-kb 200 --out ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}

# Extract the index SNP ID:
awk 'NR!=1 {print $3}' ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}.clumped >  ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}_valid_snps.txt



treatment=veh
cpg=cg12487466

echo "SNP P" > test.txt
awk -v cg=$cpg '$2 ~ cg {print $1 " " $6}' ../me-qtl_cis_result_${treatment}_fdr_005.csv >> test.txt
plink --bfile $snp_dir/Dex_genoData_SNPs --clump test.txt --clump-r2 0.2 --clump-kb 200 --clump-p1 0.05 --clump-p2 1 --out ${treatment}_${cpg}


Submitted batch job 3350908
Submitted batch job 3350909
Submitted batch job 3350910
Submitted batch job 3350911
Submitted batch job 3350912
Submitted batch job 3350913
Submitted batch job 3350914
Submitted batch job 3350915
Submitted batch job 3350916
Submitted batch job 3350917
Submitted batch job 3350918
Submitted batch job 3350919
Submitted batch job 3350920
Submitted batch job 3350921