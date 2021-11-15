cd /home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL

snp_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/snps

treatment=$1

echo "SNP P" > ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

awk 'NR != 1 {print $1 "\t" $5}' me-qtl_cis_result_${treatment}_fdr_005.csv | uniq >> ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

# sed -i '1 i\SNP\tP' ind_me-qtl_cis/me-qtl_cis_${treatment}_snps_for_ld_clump.csv

plink --bfile $snp_dir/Dex_genoData_SNPs --clump ind_me-qtl_cis/${treatment}/me-qtl_cis_${treatment}_snps_for_ld_clump.csv --clump-r2 0.2 --clump-kb 200 --out ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}

# Extract the index SNP ID:
awk 'NR!=1 {print $3}' ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}.clumped >  ind_me-qtl_cis/${treatment}/me-qtl_cis_ind_${treatment}_valid_snps.txt
