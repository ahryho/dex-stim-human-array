# pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
# setwd(pre)

library(optparse)

# CLI parsing
option_list = list(
  make_option(c("-r", "--directory"),
               action = "store",
               type = "character",
               default = NULL,
               help = "a directory where the methylation data are saved",
               metavar = "character"),
  make_option(c("-d", "--data"),
               type = "character",
               default = NULL,
               help = "an input data filename, result of lmer",
               metavar = "character"),
  make_option(c("-m", "--model"),
               type = "character",
               default = NULL,
               help = "a lmer model",
               metavar = "character"),
  make_option(c("-p", "--pval"),
                type = "double",
                default = 0.05,
                help = "p value threshold [default = %default]",
                metavar = "double"),
  make_option(c("-b", "--beta"),
                type = "double",
                default = 0.1,
                help = "delta beta threshold [default = %default]",
                metavar = "double")
)

opt <- parse_args(OptionParser(option_list=option_list))

src.data.pre    <- opt$directory
dmp.pval.adj.fn <- opt$data
mdl             <- opt$model    
delta.beta      <- opt$beta 
p.thr           <- opt$pval 

# Load libraries

source("util.R", chdir = TRUE)

pkg.list <- c("BiocManager", "tidyverse", "dplyr", "ggplot2", "Gviz", "RColorBrewer", "glue")
biocmanager.pkg.list <- c("ggman", "GenomicRanges")

cat("Loading libraries\n")
suppressPackageStartupMessages(LoadPackages(pkg.list))
suppressPackageStartupMessages(LoadPackagesBio(biocmanager.pkg.list))
cat("\nLibraries loaded\n")

# dmp.pval.adj.fn  <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_pval_adj_bcc_pcs.txt")
dmp.sign.df.fn   <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_with_beta_stat_")
dmp.sign.id.fn   <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_ids_")

report.dir       <- paste0(src.data.pre, "20_DMA/90_reports/02_dmp/")

# Load data

pval.adj.tbl          <- read.csv(dmp.pval.adj.fn, sep = "\t")
pval.adj.tbl$FC       <- as.numeric(pval.adj.tbl$FC)
pval.adj.tbl$Pval_adj <- as.numeric(pval.adj.tbl$Pval_adj)

# Get significant dmps

nr.hypermethyl <- nrow(pval.adj.tbl[pval.adj.tbl$FC > delta.beta & pval.adj.tbl$Pval_adj < p.thr,])
nr.hypomethyl  <- nrow(pval.adj.tbl[pval.adj.tbl$FC < delta.beta & pval.adj.tbl$Pval_adj < p.thr,])

pval.adj.tbl[, "threshold"] <- as.factor(abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr)
cols <- c("TRUE" = "green", "FALSE" = "grey")
pdf(file = paste0(report.dir, 
                  "Volcano_plot_model_", mdl, 
                  "beta_", delta.beta * 100, 
                  "_p_", p.thr * 100,
                  ".pdf" ))

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
dev.off()

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