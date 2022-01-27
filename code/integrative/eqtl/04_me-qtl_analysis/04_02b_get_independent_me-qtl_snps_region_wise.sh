# region-wise clumping

cd ~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/region_wise_independent_snps 

treatment=delta

wc -l me-qtl_cis_result_${treatment}_fdr_005.csv
head -n 10 me-qtl_cis_result_${treatment}_fdr_005.csv

awk 'NR != 1 {print $2}' me-qtl_cis_result_${treatment}_fdr_005.csv | sort | uniq > list_${treatment}_cpgs_for_clumping.csv 

head -n 10 list_${treatment}_cpgs_for_clumping.csv 
wc -l list_${treatment}_cpgs_for_clumping.csv 

snp_dir=~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_janine/qc/

echo "CpG_ID SNP" > me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

while read line
do
    cpg=$(echo $line | cut -f1 -d$'\t')

    echo "SNP P" > ${treatment}_snp_lst_for_ld_clump_tmp.txt
    awk -v cg=$cpg '$2 ~ cg {print $1 " " $5}' me-qtl_cis_result_${treatment}_fdr_005.csv >> ${treatment}_snp_lst_for_ld_clump_tmp.txt
    
    plink --bfile $snp_dir/Dex_genoData_SNPs --clump ${treatment}_snp_lst_for_ld_clump_tmp.txt --clump-r2 0.2 --clump-kb 200 --out ${treatment}_${cpg}
    
    awk -v cg=$cpg 'NR != 1 {if($3) print cg " " $3}' ${treatment}_${cpg}.clumped >> me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

    rm {*.log,*.clumped}
  
done < list_${treatment}_cpgs_for_clumping.csv 

# test.txt

rm ${treatment}_snp_lst_for_ld_clump_tmp.txt


wc -l ${treatment}_snp_lst_for_ld_clump_tmp.txt
head -n 10 ${treatment}_snp_lst_for_ld_clump_tmp.txt
