
pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

# Load libraries

pkg.list <- c("BiocManager", "tidyverse", "dplyr", "ggplot2", "Gviz", "RColorBrewer", "glue")
biocmanager.pkg.list <- c("ggman", "GenomicRanges")

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

dmp.pval.adj.fn  <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_pval_adj_bcc_pcs.txt")
dmp.sign.df.fn   <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_with_beta_stat_")
dmp.sign.id.fn   <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_ids_")

# Load data

pval.adj.tbl          <- read.csv(dmp.pval.adj.fn, sep = "\t")
pval.adj.tbl$FC       <- as.numeric(pval.adj.tbl$FC)
pval.adj.tbl$Pval_adj <- as.numeric(pval.adj.tbl$Pval_adj)

# Get significant dmps
mdl        <- "bcc_pcs_"    
delta.beta <- 0.1
p.thr      <- 0.05

nr.hypermethyl <- nrow(pval.adj.tbl[pval.adj.tbl$FC > delta.beta & pval.adj.tbl$Pval_adj < p.thr,])
nr.hypomethyl  <- nrow(pval.adj.tbl[pval.adj.tbl$FC < delta.beta & pval.adj.tbl$Pval_adj < p.thr,])

pval.adj.tbl[, "threshold"] <- as.factor(abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr)
cols <- c("TRUE" = "green", "FALSE" = "grey")
ggplot(pval.adj.tbl, aes(y = -log10(Pval_adj), x = FC, color = threshold)) +
  geom_point(alpha = .5, size = 1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = delta.beta, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -delta.beta, colour="#990000", linetype="dashed") +
  geom_hline(yintercept = -log10(p.thr), colour="#990000", linetype="dashed") +
  theme(legend.position = "none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste0("Volcano plot\nP-value threshold = ", p.thr, "\n", nr.hypermethyl, " hypermethyl CpGs\n", nr.hypomethyl, " hypomethyl CpGs"))

# Get significant DMPs

dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, 1:3]
dmp.sign.df  <- dmp.sign.df[order(dmp.sign.df$Pval_adj),]
dmp.sign.ids <- dmp.sign.df$Probe_Id 

write.table(dmp.sign.df,
            file = paste0(dmp.sign.df.fn, 
                          mdl, 
                          "beta_", delta.beta * 100, 
                          "_p_", p.thr * 100, 
                          ".txt"), row.names = F, quote = F, sep = "\t", col.names = T)

write.table(dmp.sign.ids,
            file = paste0(dmp.sign.fn,
                          mdl, 
                          "beta_", delta.beta * 100, 
                          "_p_", p.thr * 100, 
                          ".txt"), row.names = F, quote = F, sep = "\t", col.names = F)