
library(limma)
library(DMRcate)
library(Gviz)

library(ggman)

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])
lmer.res.out.fn     <- as.character(args[2])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

beta.stat.fn     <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_no_cov_with_statistics.txt"
# lmer.mdl.1.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_no_cov.txt"
lmer.mdl.2.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_sex.txt"
lmer.mdl.3.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_all.txt"
lmer.mdl.4.fn    <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/30_lme_models/lme_dex_all_plus_cell_counts.txt"

dmp.sign.df.fn   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_with_beta_stat.txt"
dmp.sign.fn      <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_ids.txt"
dmp.sign.anno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_annotated.txt"
dmp.sign.anno.genes.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmps_significant_annotated_genes.txt"
dmr.sign.anno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/40_significant/dmrs_significant_annotated.rds"
dmr.dmrcate.fn   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMR/dmr_by_dmrcate.rds"

# 1. Identify significant DMPs
# absolute value of delta beta (FC) > 0.15, p < 0.05

beta.statistics <- read.csv(beta.stat.fn, sep = "\t", stringsAsFactors = F)
lmer.mdl.2 <- read.table(lmer.mdl.2.fn, sep = "\t")
lmer.mdl.3 <- read.table(lmer.mdl.3.fn, sep = "\t")
lmer.mdl.4 <- read.table(lmer.mdl.4.fn, sep = "\t")

mthd <- "fdr"
pval.adj.tbl <- data.frame(Probe_Id = beta.statistics$Probe_Id, 
                           FC = beta.statistics$FC,
                           Pval_1 = p.adjust(beta.statistics$pVal_pheno, method = mthd), 
                           Pval_2 = p.adjust(lmer.mdl.2$V2, method = mthd),
                           Pval_3 = p.adjust(lmer.mdl.3$V2, method = mthd),
                           Pval_4 = p.adjust(lmer.mdl.4$V2, method = mthd)) %>% 
  rowwise() %>%
  mutate(P_adj_avg = mean(c(Pval_1, Pval_2, Pval_3, Pval_4))) %>%
  select(Probe_Id, FC, P_adj_avg)


delta.beta <- 0.15
p.thr      <- 0.05

dmp.sign.df  <- pval.adj.tbl[abs(pval.adj.tbl$FC) > delta.beta & pval.adj.tbl$P_adj_avg < p.thr, ]
dmp.sign.ids <- dmp.sign.df$Probe_Id #192 CpGs

write.table(dmp.sign.df,
            file = dmp.sign.df.fn, row.names = F, quote = F, sep = "\t", col.names = T)

write.table(dmp.sign.ids,
            file = dmp.sign.fn, row.names = F, quote = F, sep = "\t", col.names = F)

# ----------
# 
# summary(beta.statistics )
# 
# delta.beta <- 0.15
# p.thr      <- 0.05
# 
# dmp.sign.df  <- beta.statistics[abs(beta.statistics$FC) > delta.beta & beta.statistics$pVal_pheno < p.thr, ]
# dmp.sign.ids <- dmp.sign.df$Probe_Id # 216 CpGs
#   
# write.table(dmp.sign.df,
#             file = dmp.sign.df.fn, row.names = F, quote = F, sep = "\t", col.names = T)
# 
# write.table(dmp.sign.ids,
#             file = dmp.sign.fn, row.names = F, quote = F, sep = "\t", col.names = F)


# 2. Annotate the sign DMPs

dmp.anno.tbl  <- read.csv(paste0(dmr.data.dir, "dmp_annotated.csv"), sep = ";") # 740357 x 40
dmp.sign.anno <- dmp.anno.tbl[dmp.anno.tbl$Name %in% dmp.sign.ids, ] # 192 x 40

write.table(dmp.sign.anno, file = dmp.sign.anno.fn, sep = ";", row.names = F)

sign.dmp.genes <- dmp.sign.anno[dmp.sign.anno$UCSC_RefGene_Name != "",] # only those probes with assigned genes

write.table(sign.dmp.genes, file = dmp.sign.anno.genes.fn, sep = ";", row.names = F)

# 3. Plot the DMPs

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
