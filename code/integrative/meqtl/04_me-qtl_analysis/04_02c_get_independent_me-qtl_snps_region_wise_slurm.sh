#!/bin/bash
#
#SBATCH --job-name=get_ind_meSNPs
#SBATCH --output=err_out/get_ind_meSNPs_%A_%a.out
#SBATCH --error=err_out/get_ind_meSNPs_%A_%a.err
#SBATCH --mem-per-cpu=100M     # Each task uses max 2G of memory
#SBATCH --array=1-10000%100    # Submit 10'000 tasks. Run max 100 concurrently
#SBATCH --part=pe
#SBATCH --exclude=hp01,hp02,pe6


# snp_dir=~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_janine/qc/
snp_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/snps
rslt_dir=/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/region_wise_independent_snps
treatment=$2

# test.txt
# awk 'NR != 1 {print $2}' me-qtl_cis_result_${treatment}_fdr_005.csv | sort | uniq > list_${treatment}_cpgs_for_clumping.csv

# head -n 10 list_${treatment}_cpgs_for_clumping.csv
# wc -l list_${treatment}_cpgs_for_clumping.csv

# num_cpgs=$(wc -l list_${treatment}_cpgs_for_clumping.csv | awk '{print $1}')
# itrait=$(($SLURM_ARRAY_TASK_ID%$num_cpgs))

rslt_treatment_dir=${rslt_dir}/${treatment}
mkdir -p ${rslt_treatment_dir}

SCRATCH_DIRECTORY=${rslt_dir}/scratch/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# cp ${rslt_dir}/test.txt ${SCRATCH_DIRECTORY}
# cp ${rslt_dir}/me-qtl_cis_result_${treatment}_fdr_005.csv

i=$1
# i=130000
n=$(($SLURM_ARRAY_TASK_ID+${i}))
cpg=$(awk -v cpg_id=$n 'NR==cpg_id {print $1}' ${rslt_dir}/list_${treatment}_cpgs_for_clumping.csv)

echo "SNP P" > ${treatment}_snp_lst_for_ld_clump_${cpg}.txt
awk -v cg=$cpg '$2 ~ cg {print $1 " " $5}' ${rslt_dir}/me-qtl_cis_result_${treatment}_fdr_005.csv >> ${treatment}_snp_lst_for_ld_clump_${cpg}.txt

plink --bfile $snp_dir/Dex_genoData_SNPs --clump ${treatment}_snp_lst_for_ld_clump_${cpg}.txt --clump-r2 0.2 --clump-kb 200 --out ${treatment}_${cpg}

# echo "CpG_ID SNP" >  me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt
awk -v cg=$cpg 'NR != 1 {if($3) print cg " " $3}' ${treatment}_${cpg}.clumped >> ${SCRATCH_DIRECTORY}/me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt

# rm {*.log,*.clumped,${treatment}_snp_lst_for_ld_clump_${cpg}.txt}

# echo "CpG_ID SNP" > me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt
# cat ${rslt_treatment_dir}/*.txt >> ${rslt_treatment_dir}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

cp ${SCRATCH_DIRECTORY}/me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt ${rslt_treatment_dir}

# mv ${rslt_treatment_dir}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt ${rslt_dir}
# rm ${rslt_treatment_dir}/*

cd ${rslt_dir}
rm -rf ${SCRATCH_DIRECTORY}

# awk 'NR !=1 {print $1}' ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt | uniq |wc -l
# awk 'NR !=1 {print $2}' ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt | uniq |wc -l
# wc -l ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

exit 0