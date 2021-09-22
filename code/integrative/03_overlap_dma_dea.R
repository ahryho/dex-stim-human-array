source("~/bio/code/mpip/dex-stim-human-array/code/methylation/02_dma/01_dmp/funcitons.R")

library(data.table)
library(ggplot2)
library(eulerr)
library(dplyr)

# GEX

lmem.rslt.fn <- '~/bio/code/mpip/dex-stim-human-array/output/data/gene_expression/01_lmem_dea/lmem_dex_svs_1_5.csv'
lmem.rslt    <- fread(lmem.rslt.fn)

colnames(lmem.rslt) <- c("PROBE_ID", "PVAL", "CHI_SQ", 
                         "MEAN_VEH", "SD_VEH", "VAR_VEH", 
                         "MEAN_DEX", "SD_DEX", "VAR_DEX",
                         "FC", "VAR_ALL", "PVAL_PHENO", "FDR")

fdr <- 0.05
fc  <- 0.05

VolcanoPlot(lmem.rslt, fdr, fc, "GEX SVs")


deg.sign.df <- lmem.rslt[abs(FC) > fc & FDR <= fdr,]
  
fwrite(deg.sign.df, 
       "~/bio/code/mpip/dex-stim-human-array/output/data/gene_expression/01_lmem_dea/deg_fdr005_fc005.csv",
       quote = F, row.names = F)

map.ilmn.tbl.fn <- '~/bio/code/mpip/dex-stim-human-array/data/mapping/mapping_ilmn_ensg_gene.csv'
map.ilmn.gene.ensg.tbl <- fread(map.ilmn.tbl.fn)

gex.ilmn.ids <- deg.sign.df$PROBE_ID

gex.gene.ids <- map.ilmn.gene.ensg.tbl[Illumina_ID %in% gex.ilmn.ids, Gene_ID ] %>% unique()


# CpGs

map.cpg.tbl.fn  <- '~/bio/code/mpip/dex-stim-human-array/data/mapping/mapping_cpg_gene_ensg_full.csv'
map.cpg.gene.ensg.tbl  <- fread(map.cpg.tbl.fn)

dmps.sign.anno.fn <- "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dmps_fdr01_fc02_anno_full.csv"
dmps.sign.anno.df <- fread(dmps.sign.anno.fn)

meth.gene.ids <- map.cpg.gene.ensg.tbl[CpG_ID %in% dmps.sign.anno.df$PROBE_ID, Gene_ID ] %>% unique()

# Overlap

fit <- euler(list(DMP = meth.gene.ids, DEA = gex.gene.ids))
round(fit$original.value["DMP&DEA"] / length(gex.ensg.ids) * 100, 1)
round(fit$original.value["DMP&DEA"] / length(meth.ensg.ids) * 100, 1)

plot(fit,
     quantities = list(type = c("counts", "percent")),
     fill = c("red", "white"))

length(gex.ilmn.ids); length(gex.gene.ids)
length(dmps.sign.anno.df$PROBE_ID); length(meth.gene.ids)
