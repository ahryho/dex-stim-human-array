
library(limma)
library(DMRcate)

library(Gviz)
library(ggplot2)
library(RColorBrewer)

library(ggman)
library(missMethyl)

library(dplyr)
library(glue)

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])
lmer.res.out.fn     <- as.character(args[2])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
# input.parameters  <- "input_parameters.csv"
input.parameters    <- as.data.frame(input.parameters)


beta.stat.fn     <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_no_cov_with_statistics.txt"
lmer.mdl.1.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_no_cov.txt"
lmer.mdl.2.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_sex.txt"
lmer.mdl.3.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_all.txt"
lmer.mdl.4.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_all_plus_cell_counts.txt"
lmer.mdl.5.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_all_plus_cell_counts_pcs.txt"

dmp.pval.adj.fn  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/dmps_pval_adj.txt"
dmp.sign.df.fn   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_with_beta_stat.txt"
dmp.sign.fn      <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_ids.txt"
dmp.sign.anno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_annotated.txt"
dmp.sign.anno.genes.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_annotated_genes.txt"
dmr.sign.anno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmrs_significant_annotated.rds"
dmr.dmrcate.fn   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMR/dmr_by_dmrcate.rds"

# 1. Identify significant DMPs
# absolute value of delta beta (FC) > 0.15, p < 0.05

beta.statistics <- read.csv(beta.stat.fn, sep = "\t", stringsAsFactors = F)
lmer.mdl.1 <- read.table(lmer.mdl.1.fn, sep = "\t")
lmer.mdl.2 <- read.table(lmer.mdl.2.fn, sep = "\t")
lmer.mdl.3 <- read.table(lmer.mdl.3.fn, sep = "\t")
lmer.mdl.4 <- read.table(lmer.mdl.4.fn, sep = "\t")
lmer.mdl.5 <- read.table(lmer.mdl.5.fn, sep = "\t")[-1, ]

lmer.mdl <- lmer.mdl.5 
lmer.mdl$V2 <- as.numeric(lmer.mdl$V2)
lmer.mdl$V3 <- as.numeric(lmer.mdl$V3)
  
mthd <- "fdr"
pval.adj.tbl <- data.frame(Probe_Id = beta.statistics$Probe_Id, 
                           FC = beta.statistics$FC,
                           Pval_adj = p.adjust(lmer.mdl$V2, method = mthd)) 

# Compare all adj p-values

pval.adj.all.mdls <- rbind(data.frame(p_val = p.adjust(lmer.mdl.1$V2, method = mthd), mdl = "mdl1"),
                           data.frame(p_val = p.adjust(lmer.mdl.2$V2, method = mthd), mdl = "mdl2"),
                           data.frame(p_val = p.adjust(lmer.mdl.3$V2, method = mthd), mdl = "mdl3"),
                           data.frame(p_val = p.adjust(lmer.mdl.4$V2, method = mthd), mdl = "mdl4"),
                           data.frame(p_val = p.adjust(lmer.mdl.5$V2, method = mthd), mdl = "mdl5"))

pdf(file = paste0(report.dir, "PValues_adj_comparison_all_mdls_boxplot.pdf"))
ggplot(pval.adj.all.mdls, aes(x = mdl, y = p_val, color = mdl)) +
  geom_boxplot()
dev.off()

write.table(pval.adj.tbl,
            file = dmp.pval.adj.fn , row.names = F, quote = F, sep = "\t", col.names = T)

# extract significant DMPs

# estimate delta beta and p thr

# model without PCs
pval.adj.tbl <- read.csv2("/Users/anastasiia_hry/bio/workspace/dex/dex-methylation/dmps_pval_adj.txt", sep = "\t")
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

pval.adj.tbl <- read.csv2("/Users/anastasiia_hry/bio/workspace/dex/dex-methylation/dmps_pval_adj_pcs.txt", sep = "\t")
pval.adj.tbl$FC <- as.numeric(pval.adj.tbl$FC)
pval.adj.tbl$Pval_adj <- as.numeric(pval.adj.tbl$Pval_adj)

# annotate dmps
# dmp.df <- pval.adj.tbl
# 
# anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmp.df$Probe_Id, c(1:4,12:19,22:43)])
# 
# dmp.anno.df <- right_join(anno.epic.sub, dmp.df, by = c("Name" = "Probe_Id"))  # 740,357 x 36
# dmp.anno.df <- dmp.anno.df[order(dmp.sign.anno$Pval_adj),]
# 
# dir.prefix  <- "/Users/anastasiia_hry/bio/datasets/methylation/20_DMA/"
# dmp.anno.fn <- paste0(dir.prefix, "02_dmp/dmp_bcc_pcs_anno.csv")
# 
# write.table(dmp.anno.df, file = dmp.anno.fn, sep = ";", row.names = F)

############

pval.list <- c(0.01, 0.02, 0.03, 0.04, 0.05)
beta.list <- c(0, 0.01, 0.05, 0.1, 0.15, 0.17)

# pvalues
nr.dmps.list <- c()
for (p.thr in pval.list){
  # delta.beta <- 0.0
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


delta.beta <- 0.0
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

nr.hypermethyl <- nrow(pval.adj.tbl[pval.adj.tbl$FC > delta.beta & pval.adj.tbl$Pval_adj < p.thr,])
nr.hypomethyl  <- nrow(pval.adj.tbl[pval.adj.tbl$FC < delta.beta & pval.adj.tbl$Pval_adj < p.thr,])

# get significant DMPs

dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$Pval_adj < p.thr, 1:3]
dmp.sign.df  <- dmp.sign.df[order(dmp.sign.df$Pval_adj),]
dmp.sign.ids <- dmp.sign.df$Probe_Id 

write.table(dmp.sign.df,
            file = dmp.sign.df.fn, row.names = F, quote = F, sep = "\t", col.names = T)

write.table(dmp.sign.ids,
            file = dmp.sign.fn, row.names = F, quote = F, sep = "\t", col.names = F)


# 2. Annotate significant DMPs 

dmr.data.dir <- "/Users/anastasiia_hry/bio/workspace/dex/dex-methylation/"
dmp.sign.anno.fn <- paste0(dmr.data.dir, "dmps_significant_annotated.txt")

anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmp.sign.ids, c(1:4,12:19,22:43)])

dmp.sign.anno <- right_join(anno.epic.sub, dmp.sign.df, by = c("Name" = "Probe_Id"))  # 3560 x 36
dmp.sign.anno <- dmp.sign.anno[order(dmp.sign.anno$Pval_adj),]

write.table(dmp.sign.anno, file = dmp.sign.anno.fn, sep = ";", row.names = F)

# sign.dmp.genes <- dmp.sign.anno[dmp.sign.anno$UCSC_RefGene_Name != "",] # only those probes with assigned genes
# write.table(sign.dmp.genes, file = dmp.sign.anno.genes.fn, sep = ";", row.names = F)

# dmp.sign.anno <- read.csv2(dmp.sign.anno.fn)

# Plot the top few differentially methylated CpGs to check whether the results make sense

pdf(file = paste0(report.dir, "04_b_Sign_DMC_top_10.pdf"))
par(mfrow = c(2,2))
sapply(dmp.sign.df$Probe_Id[1:10], function(cpg){
  plotCpg(beta.mtrx, cpg = cpg, pheno = pheno$Group, ylab = "Beta values")
})
dev.off()

# 3. Flat annotation

flatAnn <- missMethyl:::.getFlatAnnotation(anno.epic.sub, array.type = "EPIC")

# Associate CpGs to genes (ENTREZ ID) using the Illumina annotation information.
GetMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

dat <- data.frame(table(table(flatAnn$entrezid)))
numCpgsPerGene <- as.vector(table(flatAnn$entrezid))
med <- median(numCpgsPerGene)
mod <- GetMode(numCpgsPerGene)

# Frequencies per gene

ggplot(dat, aes(x = Var1, y = Freq, group = 1)) +
  geom_line(linetype = "dotted", size = 1) +
  geom_point( color="black", size = 2) +
  geom_vline(xintercept = med, linetype = "dashed", color = "red") +
  annotate("text", x = med + 0.1, label = glue("median = {med}"), y = 720, colour = "red", 
           size = 3, hjust="left") +
  geom_vline(xintercept = mod, linetype = "dashed", color = "blue") +
  annotate("text", x = mod + 0.1, label = glue("mode = {mod}"), y = 780, colour = "blue", 
           size = 3, hjust="left") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") +
  xlab("No. CpGs per gene") +
  ylab("Frequency") +
  geom_text(aes(label = Freq), vjust = -0.5, nudge_y = 0.) 

# Frequencies per gene per GO group

dat <- data.frame(table(flatAnn$group))
dat <- rbind(dat, data.frame(Var1 = "Intergenic", Freq = nrow(dmp.sign.df) - nrow(flatAnn)))
dat$Freq <- round(dat$Freq / nrow(dmp.sign.anno) * 100, 2)
ggplot(dat, aes(x = Var1, y = Freq, group = 1)) +
  geom_bar(stat="identity") +
  xlab("No. CpGs per gene per GO category") +
  ylab("Frequency") +
  geom_text(aes(label = Freq), vjust = -0.5, nudge_y = 0.) 

# 4. Significant associations between DNA methylation and Group

flat.ann.pval <- left_join(flatAnn, dmp.sign.df, by = c("cpg" = "Probe_Id"))
flat.ann.pval <- flat.ann.pval[order(flat.ann.pval$Pval_adj), ]

# 4. Check the SNPs (significant associations with CpGs) overlap with Janine 3662 SNPs 

grsnp.fn <- "/Users/anastasiia_hry/bio/workspace/dex/dex-genotype/02_GR_SNPs_hg19.txt"
grsnp.df <- read.csv(grsnp.fn, sep = "\t")

dmp.snp.df <- dmp.sign.anno[!is.na(dmp.sign.anno$Probe_rs),] #532
dmp.snp.df[dmp.snp.df$Probe_rs %in% grsnp.df$ID,c(1, 2, 4, 5, 6, 35, 36)]
View(dmp.snp.df[,c(1, 2, 4, 5, 6, 35, 36)])

# chr       pos       Name   Probe_rs Probe_maf     FC     Pval_adj
# chr13 114828455 cg16962442 rs61973909  0.099634 -0.101 0.0006270184
# chr7  43676945 cg06012582 rs17439918  0.100073  0.086 0.0295883821

library(ggVennDiagram)
ggVennDiagram(list(DMPs = dmp.snp.df$Probe_rs, GR_SNPs = grsnp.df$ID)) + theme(legend.position = "none")

# 7. Plot the PCA DMPs

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

beta.mtrx <- readRDS(beta.mtrx.fn)
mval.mtrx <- readRDS(mval.mtrx.fn)

pheno.methyl.pcs.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMR/pheno_with_pcs.csv"
pheno     <- read.csv2(pheno.methyl.pcs.fn)

pal <- brewer.pal(8,"Set2")
pdf(file = paste0(report.dir, "PCA_mval_sign_cpgs_Group.pdf"))
plotMDS(mval.mtrx[dmp.sign.ids,], gene.selection = "common", 
        col = pal[factor(pheno$Group)], labels = pheno$Group, 
        main = " ", cex = 0.5)
# legend("bottom", legend = levels(factor(pheno$Group)), text.col = pal)
dev.off()


# 5. GO

sig.cpgs <- dmp.sign.anno$Name
all.cpgs <- dmp.anno.tbl$Name

gst <- gometh(sig.cpg = sig.cpgs, all.cpg = all.cpgs, plot.bias = TRUE)
topGSA(gst, number=10)


pdf(file = paste0(report.dir, "PCA_mval_sign_cpgs_Sex.pdf"))
plotMDS(mval.mtrx[dmp.sign.ids,], gene.selection = "common", 
        col = pal[factor(pheno$Sex)], labels = pheno$Sex, 
        main = " ", cex = 0.5)
# legend("bottom", legend = levels(factor(pheno$Group)), text.col = pal)
dev.off()















dmp.anno.tbl  <- read.csv2("/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMR/dmp_annotated.csv", sep = ";")
dmp.sign.anno <- read.csv2("/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_annotated.txt", sep = ";")

# dmp.sign.anno <- read.csv2("/Users/anastasiia_hry/bio/workspace/dex-methylation/dmps_significant_annotated.txt", sep = ";")
# dmp.anno.tbl  <- read.csv2("/Users/anastasiia_hry/bio/workspace/dex-methylation/dmp_annotated.csv", sep = ";")

df   <- dmp.anno.tbl
CHR  <- as.numeric(substr(df$chr, 4, length(df$chr)))
BP   <- as.numeric(df$pos)

P    <- as.numeric(df$P.Value)
P[is.na(P)] <- 1

SNP  <- df$Probe_rs

GENE <- df$GencodeCompV12_NAME
GENE <- gsub( ".*;", "", GENE)

df <- cbind.data.frame(CHR,BP,P, SNP, GENE)

# ---

df.sign <- dmp.sign.anno
CHR     <- as.numeric(substr(df.sign$chr, 4, length(df.sign$chr)))
BP      <- as.numeric(df.sign$pos)

P    <- as.numeric(df.sign$P.Value)
P[is.na(P)] <- 1

SNP  <- df.sign$Probe_rs

GENE <- df.sign$GencodeCompV12_NAME
GENE <- gsub( ".*;", "", GENE)

df.sign <- cbind.data.frame(CHR,BP,P, SNP, GENE)

# ---

th <- max(as.numeric(dmp.sign.anno$P.Value))
snp.of.interest <- na.omit(dmp.sign.anno$Probe_rs)

# Label genes
genes.id.top10 <- dmp.sign.anno$GencodeCompV12_NAME
genes.id.top10 <- na.omit(dmp.sign.anno$GencodeCompV12_NAME)
genes.id.top10 <- genes.id.top10[genes.id.top10 != ""]
genes.id.top10 <- gsub( ".*;", "", genes.id.top10 )


# manhattan(df, cex = 1, colors = c("blue","green"), suggestiveline = -log10(th), highlight = snp.of.interest, cex.axis = 1, main = '')

pdf(file = "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/Manhattan_plot_dmps.pdf")
clumps <- ggmanClumps(df.sign, index.snp.column = "SNP", clumps.column = "GENE", label.column = "GENE") 
dmps.pl.man <- ggman(df[1:1000, ], clumps = clumps, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P", sigLine = -log10(th), relative.positions = T)

dmps.pl.man <- ggman(df, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P", sigLine = -log10(th), relative.positions = T)
ggmanLabel(dmps.pl.man, labelDfm = df.sign[1:10, ], snp = "SNP", label = "GENE", type = 'text')
dev.off()

# 3. Significant DMR

dmrs <- readRDS(dmr.dmrcate.fn) # 54762 ranges

results.ranges <- dmrs

gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

ann450kSub <- dmp.sign.anno 
ann450kOrd <- ann450kSub[order(ann450kSub$chr,ann450kSub$pos),]
bValsOrd <- beta.mtrx[match(ann450kOrd$Name,rownames(beta.mtrx)),]

cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bValsOrd)

# extract data on CpGs in DMR

cpgDataDMR <- subsetByOverlaps(cpgData, dmrs[dmrIndex])

# 4. Plot








pheno      <- read.csv2(pheno.fn)

Group      <- factor(pheno$Group)
Individual <- factor(pheno$Individual)

design           <- model.matrix(~ 0 + Group + Individual, data = pheno)
colnames(design) <- c(levels(Group), levels(Individual)[-1])

contMatrix <- makeContrasts(dex-veh, levels = design)
contMatrix

# 3. Annotation for M values 

mval.mtrx   <- readRDS(mval.mtrx.fn)
mval.subset <- mval.mtrx[rownames(mval.mtrx) %in% dmp.sign.ids, ] # 216 x 399

# anno <- readRDS(paste(dmr.data.dir, "anno_for_dmr.rds"))
anno <- cpg.annotate(object = mval.subset, datatype = "array", what = "M", 
                     analysis.type = "differential", design = design, 
                     contrasts = TRUE, cont.matrix = contMatrix, 
                     coef = "dex - veh", arraytype = "EPIC")

# Once we have the relevant statistics for the individual CpGs, 
# we can then use the dmrcate function to combine them to identify differentially methylated regions. 
# The main output table dmr.tbl$results contains all of the regions found, 
# along with their genomic annotations and p-values. 

dmr.tbl        <- dmrcate(anno, lambda = 1500, C = 2, pcutoff = 0.95) # sigma = lambda / C
results.ranges <- extractRanges(dmr.tbl)
results.ranges
saveRDS(results.ranges, dmr.sign.anno.fn)
