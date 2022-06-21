
library(data.table)
library(dplyr)

## ---------
# Main function
## ---------

GetMethylChange <- function(meqtl.select, snp.df, methyl.df, treatment = "dex"){
  
  snp.id <- meqtl.select[1] #"rs79375669" # "rs10460588"
  cpg.id <- meqtl.select[2] #"cg02184558" # 
  
  snp.layer <- snp.df[SNP == snp.id]
  methyl.layer <- methyl.df[CpG_ID == cpg.id]
  
  meqtl <- data.table(t(rbind(snp.layer[, -1], methyl.layer[, -1])))
  colnames(meqtl) <- c("SNP", "CpG")
  
  cnt.df <- meqtl[, .(mean = mean(CpG)), by = SNP]
  
  if(nrow(cnt.df) == 3){
    cnt.df <- cnt.df[SNP == 0, col := "ref_mean"]
    cnt.df <- cnt.df[SNP == 1, col := "heterozyg_mean"]
    cnt.df <- cnt.df[SNP == 2, col := "non_ref_mean"]
    
    res <- data.table(t(cnt.df$mean))
    colnames(res) <- cnt.df$col 
    res[, methyl_change := non_ref_mean - ref_mean]
    
    res[, SNP := snp.layer[, 1]]
    res[, CpG_ID := methyl.layer[, 1]]
    res[, treatment := treatment]
  }
}

## ---------
# Set up parameters
## ---------

# eqtm.in.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
eqtm.in.pre  <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"

# out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"
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

methyl.change.dex.df <- do.call("rbind", apply(meqtl.ind.df, 1, GetMethylChange, snp.df, methyl.beta.df, treatment))

write.table(methyl.change.dex.df, file = paste0(out.dir.pre, "me-qtl_cis_indp_methyl_change_", treatment, "_fdr_005.csv"), row.names = F, quote = F, sep = "\t", col.names = T, append = F)

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

methyl.change.veh.df <- do.call("rbind", apply(meqtl.ind.df, 1, GetMethylChange, snp.df, methyl.beta.df, treatment))

write.table(methyl.change.veh.df,
            file = paste0(out.dir.pre, "me-qtl_cis_indp_methyl_change_", treatment, "_fdr_005.csv"), row.names = F, quote = F, sep = "\t", col.names = T, append = F)