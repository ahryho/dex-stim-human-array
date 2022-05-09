library(data.table)
library(dplyr)

setwd("~/bio/code/mpip/dex-stim-human-array/")

source("code/util.R")

# Load publicly available data from Josine L.Min at el. 2021

# mqtl_data <- fromJSON("http://api.godmc.org.uk/v0.1/assoc_meta/cpg/cg17242362", simplifyVector = FALSE)
 
public.meqtl.df <- fread("data/public_data/GoDMC_mQTL/assoc_meta_all.csv")
public.meqtl.cis.sub.df <- public.meqtl.df[cistrans == T, .(cpg, snp, pval)] # 71'616'458
public.meqtl.cis.sub.df[["meQTL_ID"]] <- paste(public.meqtl.cis.sub.df$cpg, public.meqtl.cis.sub.df$snp, sep = "-")

fwrite(public.meqtl.cis.sub.df, 
       "data/public_data/GoDMC_mQTL/assoc_meta_cis_meqtl.csv",
       quote = F, row.names = F, sep = ";")

public.meqtl.cis.sub.df <- fread("data/public_data/GoDMC_mQTL/assoc_meta_cis_meqtl.csv")

# Genotype data

snp.loc.df <- fread("data/integrative/matrixEQTL/snp_locations.csv")

head(snp.loc.df)
dim(snp.loc.df)

snp.loc.df[["SNP_ID"]] <- paste0("chr", snp.loc.df$chr, ":", snp.loc.df$pos, ":", "SNP")

# Baseline meQTLs
col.names    <- c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "FDR")

meqtl.veh.df <- fread("output/data/integrative/matrixEQTL/meqtls/me-eqtl_cis_results_06122021/me-qtl_cis_result_veh_fdr_005.csv", col.names = col.names) 

meqtl.veh.nom.df <- fread("output/data/integrative/matrixEQTL/meqtls/me-eqtl_cis_results_06122021/me-qtl_cis_result_veh.csv", col.names = col.names) 

# Join SNP_ID tp meQTL baseline

meqtl.veh.df  <- left_join(meqtl.veh.df, snp.loc.df)
meqtl.veh.lst <- paste(meqtl.veh.df$CpG_ID, meqtl.veh.df$SNP_ID, sep = "-") # 6'944'026

# Overlap at meQTL level
# 
# Significant @ FDR <= 0.05
intersect.meqtls <- intersect(meqtl.veh.lst, public.meqtl.cis.sub.df$meQTL_ID) # 2'851'678

fwrite(x = as.data.frame(public.meqtl.cis.sub.df$meQTL_ID), 
       "data/public_data/GoDMC_mQTL/assoc_meta_cis_meqtl_lst.csv",
       quote = F, row.names = F, sep = ";")

length(intersect.meqtls) / length(meqtl.veh.lst) # 41.1 %


# Overlap at meCpG level
# 
veh.mecpgs <- meqtl.veh.df$CpG_ID %>% unique() # 137'537
pub.mecpgs <- public.meqtl.cis.sub.df$cpg %>% unique() # 232'477

intersect.mecpgs <- intersect(veh.mecpgs, pub.mecpgs) #53'891

length(intersect.mecpgs) / length(veh.mecpgs) # 39.18%
length(intersect.mecpgs) / length(pub.mecpgs) # 23.18 %

# All 740K CpGs
# 
methyl.beta.mtrx <- LoadMethylBeta("veh")
all.cpgs         <- methyl.beta.mtrx$CpG_ID

intersect.all.cpgs <- intersect(all.cpgs, pub.mecpgs) # 203'696

length(intersect.mecpgs) / length(intersect.all.cpgs)  # 26.45%
