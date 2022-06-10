# prepare dnam beta matrix: split the entire mtric into 22 matrices: for each autosomal chromosome
# 
library(data.table)
library(dplyr) 

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/")

treatment <- "dex"

# Beta mtrx
# 
methyl.mtrx.veh <- fread(paste0("data/integrative/matrixEQTL/methyl_beta_mtrx_", treatment, ".csv"))

# Genome coordanates for CpGs
#
cpg.loc         <- fread("data/integrative/matrixEQTL/cpg_locations.csv") 

# Split by chromosome:
# 
invisible(
  lapply(unique(cpg.loc$chr), function(chr.i){
  cpg.lst <- cpg.loc[chr == chr.i, CpG_ID]
  beta.mtrx.chr <- methyl.mtrx.veh[CpG_ID %in% cpg.lst,]
  fwrite(beta.mtrx.chr,
         paste0("data/methylation/beta_mtrx_by_chrom/methyl_beta_mtrx_", treatment, "_chr_", chr.i, ".csv"),
         quote = F, row.names = F, sep = "\t")
  print(paste0("Chromosome ", chr.i, " done"))
}))
