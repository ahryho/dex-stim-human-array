# Check for overlaps between delta meQTL SNps and GR-meSNPs from Moore (LD clumped and significant at FDR = 5%)
# 
# 
setwd("~/bio/code/mpip/dex-stim-human-array/")
source("code/integrative/util.R")

library(data.table)
library(dplyr)
library(eulerr)

public.eqtl.df <- fread("data/public_data/moore_cisresults_dexStim_emp2.txt")
meqtl.delta.df <- fread("output/data/integrative/matrixEQTL/meqtls/me-qtl_cis_result_delta_fdr_005.csv") 

summary(public.eqtl.df)

public.esnps <- public.eqtl.df$SNP %>% unique()
delta.mesnp  <- meqtl.delta.df$SNP %>% unique()

euler.fit <- euler(list(`GR-eQTL-SNPs` = public.esnps,
                        `GR-meQTL-SNPs` = delta.mesnp))

plot(euler.fit,
     fills = list(fill = c("grey", "skyblue"), alpha = 0.7),
     edges = list(col = "white", lex = 2),
     quantities = list(type = c("percent", "counts"), cex = 1))

