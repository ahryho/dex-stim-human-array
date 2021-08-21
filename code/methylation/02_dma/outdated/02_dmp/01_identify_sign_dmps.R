
pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

# Load libraries

pkg.list <- c("BiocManager", "tidyverse", "dplyr", "ggplot2", "Gviz", "RColorBrewer", "glue")
biocmanager.pkg.list <- c("ggman", "missMethyl", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "GenomicRanges")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# args                <- commandArgs(T)
# input.parameters.fn <- as.character(args[1])
# lmer.res.out.fn     <- as.character(args[2])
# input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )

input.parameters  <- read.csv2(paste0(pre, "input_parameters.csv"))
input.parameters  <- as.data.frame(input.parameters)

# Set up parameters

# src.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/" 
src.data.pre <- "/Users/anastasiia_hry/bio/datasets/methylation/" 

beta.stat.fn     <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_no_cov_with_statistics.txt")
lmer.mdl.1.fn    <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_no_cov.txt")
lmer.mdl.2.fn    <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_sex.txt")
lmer.mdl.3.fn    <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_all.txt")
lmer.mdl.4.fn    <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_all_plus_cell_counts.txt")
lmer.mdl.5.fn    <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_all_plus_cell_counts_pcs.txt")

dmp.pval.adj.bcc.pcs.fn  <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_pval_adj_bcc_pcs.txt")
dmp.pval.adj.bcc.fn      <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_pval_adj_bcc.txt")

dmp.sign.anno.fn <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_annotated.txt")
dmp.sign.anno.genes.fn <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_annotated_genes.txt")
dmr.sign.anno.fn <- paste0(src.data.pre, "20_DMA/02_dmp/dmrs_significant_annotated.rds")

# 1. Identify significant DMPs
# absolute value of delta beta (FC) > 0.15, p < 0.05

beta.statistics <- read.csv(beta.stat.fn, sep = "\t", stringsAsFactors = F)
lmer.mdl.1 <- read.table(lmer.mdl.1.fn, sep = "\t")
lmer.mdl.2 <- read.table(lmer.mdl.2.fn, sep = "\t")
lmer.mdl.3 <- read.table(lmer.mdl.3.fn, sep = "\t")
lmer.mdl.4 <- read.table(lmer.mdl.4.fn, sep = "\t")
lmer.mdl.5 <- read.table(lmer.mdl.5.fn, sep = "\t")[-1, ]

# Compare all adj p-values
mthd <- "fdr"
pval.adj.all.mdls <- rbind(data.frame(p_val = p.adjust(lmer.mdl.1$V2, method = mthd), mdl = "mdl1"),
                           data.frame(p_val = p.adjust(lmer.mdl.2$V2, method = mthd), mdl = "mdl2"),
                           data.frame(p_val = p.adjust(lmer.mdl.3$V2, method = mthd), mdl = "mdl3"),
                           data.frame(p_val = p.adjust(lmer.mdl.4$V2, method = mthd), mdl = "mdl4"),
                           data.frame(p_val = p.adjust(lmer.mdl.5$V2, method = mthd), mdl = "mdl5"))

pdf(file = paste0(report.dir, "PValues_adj_comparison_all_mdls_boxplot.pdf"))
ggplot(pval.adj.all.mdls, aes(x = mdl, y = p_val, color = mdl)) +
  geom_boxplot()
dev.off()

# Save 

lmer.mdl <- lmer.mdl.5 # bcc + pcs
lmer.mdl$V2 <- as.numeric(lmer.mdl$V2)
lmer.mdl$V3 <- as.numeric(lmer.mdl$V3)

mthd <- "fdr"
pval.adj.tbl <- data.frame(Probe_Id = beta.statistics$Probe_Id, 
                           FC = beta.statistics$FC,
                           Pval_adj = p.adjust(lmer.mdl$V2, method = mthd)) 
write.table(pval.adj.tbl,
            file = dmp.pval.adj.bcc.pcs.fn , row.names = F, quote = F, sep = "\t", col.names = T)

# extract significant DMPs

# estimate delta beta and p thr

# model without PCs
pval.adj.tbl <- read.csv2(dmp.pval.adj.bcc.fn, sep = "\t")
pval.adj.tbl$FC <- as.numeric(pval.adj.tbl$FC)
pval.adj.tbl$Pval_adj <- as.numeric(pval.adj.tbl$Pval_adj)

pval.list <- c(0.01, 0.02, 0.03, 0.04, 0.05)
beta.list <- c(0, 0.01, 0.05, 0.1, 0.15, 0.17)

# pvalues
nr.dmps.list <- c()
for (p.thr in pval.list){
  delta.beta <- 0.0
  dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, ]
  nr.dmps.list <- append(nr.dmps.list, dim(dmp.sign.df)[1])
}

pval.est <- as.data.frame(cbind(Pval = pval.list, NrDMPs = nr.dmps.list, PrDMPs = round(nr.dmps.list / nrow(pval.adj.tbl) * 100, 2), mdl = "noPCs"))

# delta beta
nr.dmps.list <- c()
for (delta.beta in beta.list){
  p.thr <- 0.05
  dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, ]
  nr.dmps.list <- append(nr.dmps.list, dim(dmp.sign.df)[1])
}

delta.beta.est <- as.data.frame(cbind(Delta_beta = beta.list, NrDMPs = nr.dmps.list, PrDMPs = round(nr.dmps.list / nrow(pval.adj.tbl) * 100, 2), mdl = "noPCs"))

# model with PCs

pval.adj.tbl <- read.csv2(dmp.pval.adj.bcc.pcs.fn, sep = "\t")
pval.adj.tbl$FC <- as.numeric(pval.adj.tbl$FC)
pval.adj.tbl$Pval_adj <- as.numeric(pval.adj.tbl$Pval_adj)

# Tune p and beta threshold

pval.list <- c(0.01, 0.02, 0.03, 0.04, 0.05)
beta.list <- c(0, 0.01, 0.05, 0.1, 0.15, 0.17)

# pvalues
nr.dmps.list <- c()
for (p.thr in pval.list){
  delta.beta <- 0.0
  dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, ]
  nr.dmps.list <- append(nr.dmps.list, dim(dmp.sign.df)[1])
}

pval.est.pcs <- as.data.frame(cbind(Pval = pval.list, NrDMPs = nr.dmps.list, PrDMPs = round(nr.dmps.list / nrow(pval.adj.tbl) * 100, 2), mdl = "PCs"))
pval.est <- rbind(pval.est, pval.est.pcs)
pval.est[, 1:3] <- apply(pval.est[, 1:3], 2, as.numeric)

# delta beta
nr.dmps.list <- c()
for (delta.beta in beta.list){
  p.thr <- 0.05
  dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, ]
  nr.dmps.list <- append(nr.dmps.list, dim(dmp.sign.df)[1])
}

delta.beta.est.pcs <- as.data.frame(cbind(Delta_beta = beta.list, NrDMPs = nr.dmps.list, PrDMPs = round(nr.dmps.list / nrow(pval.adj.tbl) * 100, 2), mdl = "PCs"))
delta.beta.est <- rbind(delta.beta.est, delta.beta.est.pcs)
delta.beta.est[, 1:3] <- apply(delta.beta.est[, 1:3], 2, as.numeric)

# plot
ggplot(pval.est, aes(x = Pval, y = NrDMPs)) + 
  geom_line(aes(colour = mdl), linetype="dotted", size=2) +
  geom_point(aes(colour = mdl), size = 3) +
  geom_label(aes(label = NrDMPs), vjust = 0, nudge_y = 0) +
  theme(legend.position = "bottom")
#  geom_text(aes(label = PrDMPs), nudge_y = - 180)

ggplot(delta.beta.est, aes(x = Delta_beta, y = NrDMPs)) + 
  geom_line(aes(colour = mdl), linetype="dotted", size=2) +
  geom_point(aes(colour = mdl), size = 3) + 
  geom_label(aes(label = NrDMPs), hjust = 0.2, nudge_y = 0) +
  theme(legend.position = "bottom")
# geom_text(aes(label = PrDMPs), nudge_y = -90)


