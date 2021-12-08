
library(data.table)
library(dplyr)

## ---------
# Main function
## ---------

GetMethylChange <- function(snp, meqtl.df, methyl.beta.df, treatment = "dex"){
  snp.sub.df       <- snp.df[SNP %in% snp]
  meqtl.snp.cpg.id <- meqtl.df[SNP %in% snp.sub.df$SNP]
  
  sample.id.ref   <- colnames(snp.sub.df)[snp.sub.df == 0]
  homozyg.ref <- methyl.beta.df[rownames(methyl.beta.df) %in% meqtl.snp.cpg.id$CpG_ID, .(ref_mean = rowMeans(.SD)), .SDcols = sample.id.ref]
  
  sample.id.non.ref       <- colnames(snp.sub.df)[snp.sub.df == 2]
  homozyg.non.ref <- methyl.beta.df[rownames(methyl.beta.df) %in% meqtl.snp.cpg.id$CpG_ID, .(non_ref_mean = rowMeans(.SD)), .SDcols = sample.id.non.ref]
  
  sample.id.heter <- colnames(snp.sub.df)[snp.sub.df == 1]
  heterozyg <- methyl.beta.df[rownames(methyl.beta.df) %in% meqtl.snp.cpg.id$CpG_ID, .(heterozyg_mean = rowMeans(.SD)), .SDcols = sample.id.heter]
  
  res <- cbind(meqtl.snp.cpg.id, heterozyg, homozyg.ref, homozyg.non.ref)
  
  res[is.na(ref_mean), ref_mean := 0 ]
  res[is.na(non_ref_mean), non_ref_mean := 0 ]
  
  res[, methyl_change := ref_mean - non_ref_mean]
  res[, ref_mean_darina := (2 * ref_mean + heterozyg_mean) / (2 * (length(sample.id.ref) + length(sample.id.heter)))]
  res[, non_ref_mean_darina := (2 * non_ref_mean + heterozyg_mean) / (2 * (length(sample.id.non.ref) + length(sample.id.heter)))]
  res[, methyl_change_darina := ref_mean_darina - non_ref_mean_darina]
  res[, treatment := treatment]
}

## ---------
# Set up parameters
## ---------

eqtm.in.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
eqtm.in.pre  <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"
out.dir.pre  <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

snp.df.fn <- paste0(eqtm.in.pre, "snp_mtrx.csv")
snp.df    <- fread(snp.df.fn)

## ---------
# Dexamethasone
## ---------

treatment <- "dex"

snp.lst.fn   <- paste0(out.dir.pre, "ind_me-qtl_cis/", treatment, "/me-qtl_cis_ind_", treatment, "_valid_snps.txt")
snp.lst      <- fread(snp.lst.fn, col.names = "SNP")

meqtl.fn  <- paste0(out.dir.pre, "me-qtl_cis_result_", treatment, "_fdr_005.csv")
col.names <-  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")
meqtl.df  <- fread(meqtl.fn, col.names = col.names)
meqtl.ind.df  <- meqtl.df[SNP %in% snp.lst$SNP]

methyl.beta.df.fn <- paste0(eqtm.in.pre, "methyl_beta_mtrx_", treatment, ".csv")
methyl.beta.df    <- fread(methyl.beta.df.fn) 
methyl.beta.df    <- methyl.beta.df[CpG_ID %in% meqtl.ind.df$CpG_ID]
rownames(methyl.beta.df) <- methyl.beta.df$CpG_ID

methyl.change.dex.df <- do.call("rbind", apply(snp.lst, 1, GetMethylChange, meqtl.ind.df, methyl.beta.df, treatment))
# methyl.change.dex.df[is.na(methyl.change.dex.df$ref_mean), "ref_mean" ] <- 0

write.table(methyl.change.dex.df, file = paste0(out.dir.pre, "me-qtl_cis_ind_methyl_change_", treatment, "_fdr_005.csv"), row.names = F, quote = F, sep = "\t", col.names = T, append = F)

x <- do.call("rbind", apply(snp.lst[1], 1, GetMethylChange, meqtl.ind.df, methyl.beta.df, treatment))

## ---------
# Baseline
## ---------

treatment <- "veh"

snp.lst.fn   <- paste0(out.dir.pre, "ind_me-qtl_cis/", treatment, "/me-qtl_cis_ind_", treatment, "_valid_snps.txt")
snp.lst      <- fread(snp.lst.fn, col.names = "SNP")

meqtl.fn  <- paste0(out.dir.pre, "me-qtl_cis_result_", treatment, "_fdr_005.csv")
col.names <-  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")
meqtl.df  <- fread(meqtl.fn, col.names = col.names)

meqtl.ind.df  <- meqtl.df[SNP %in% snp.lst$SNP]

methyl.beta.df.fn <- paste0(eqtm.in.pre, "methyl_beta_mtrx_", treatment, ".csv")
methyl.beta.df    <- fread(methyl.beta.df.fn) 
methyl.beta.df    <- methyl.beta.df[CpG_ID %in% meqtl.ind.df$CpG_ID]
rownames(methyl.beta.df) <- methyl.beta.df$CpG_ID

methyl.change.veh.df <- do.call("rbind", apply(snp.lst, 1, GetMethylChange, meqtl.ind.df, methyl.beta.df, treatment))

write.table(methyl.change.veh.df,
            file = paste0(out.dir.pre, "me-qtl_cis_ind_methyl_change_", treatment, "_fdr_005.csv"), row.names = F, quote = F, sep = "\t", col.names = T, append = F)