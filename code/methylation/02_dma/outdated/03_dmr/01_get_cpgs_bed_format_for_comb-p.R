
pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

# Load libraries

pkg.list <- c("BiocManager", "tidyverse", "dplyr")
LoadPackages(pkg.list)

# Input params

# src.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/" 
src.data.pre <- "/Users/anastasiia_hry/bio/datasets/methylation/" 
dmps.anno.fn <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_annotated.csv") 
lmem.rslt.fn <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_svs.txt")

# out.data.dir <-paste0(dir.prefix, "03_dmr/")

# Load data
lmem.rslt    <- fread(lmem.rslt.fn, sep = "\t") 
dmps.anno.df <- read.table(dmps.anno.fn, header = T, sep = "\t")

# Calculate FDR

lmem.rslt[["FDR_BH"]]  <- p.adjust(lmem.rslt[["PVAL"]], method = "fdr")
lmem.rslt[["FDR_BON"]] <- p.adjust(lmem.rslt[["PVAL"]], method = "bonferroni")

# Transform the DMP df into bed format : chr | start | end | p| probe

dmps.anno.df <- right_join(dmps.anno.df, lmem.rslt[PVAL <= 0.5,], by = c("Name" = "PROBE_ID")) 
dmp.bed.df   <- dmps.anno.df %>% mutate(chrom = chr,#chr = as.numeric(substr(chr, 4, length(chr))),
                                      start = pos, 
                                      end = start + 1,
                                      p = FDR_BON, 
                                      probe = Name)  %>% 
  dplyr::select(chrom, start, end, p, probe)

dmp.bed.df$start <- as.numeric(dmp.bed.df$start)
dmp.bed.df$end <- as.numeric(dmp.bed.df$end)
dmp.bed.df$p <- as.numeric(dmp.bed.df$p)
dmp.bed.df <- dmp.bed.df[order(dmp.bed.df$chrom, dmp.bed.df$start),]

mdl <- "svs" 

write.table(dmp.bed.df, 
            paste0(src.data.pre, "02_dmp/dmps_", mdl, "_for_combp.bed"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

# p.thr <- 0.1
# 
# write.table(dmp.bed.df[dmp.bed.df$p < p.thr,], 
#             paste0(src.data.pre , "20_DMA/02_dmp/dmp_", mdl, "p_", p.thr * 100, ".bed"), 
#             col.names = T, row.names = F, sep = "\t", quote = F)
