
source("~/bio/code/mpip/dex-stim-human-array/code/util.R", chdir = TRUE)

pkg.list <- c("tidyverse", "dplyr", "data.table")
biocmanager.pkg.list <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# Set up params

dmps.anno.fn      <- '~/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dex_cpgs_annotated.csv'
dmps.sign.anno.fn <- "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dmps_fdr01_fc02_anno_full.csv"

# Load data
dmps.anno.df <- fread(dmps.anno.fn, sep = "\t")
dmps.sign.anno.df <- fread(dmps.sign.anno.fn)

colnames(dmps.anno.df)[4] <- "PROBE_ID"

GetDistances <- function(dmps.anno.df, ensembl){
  
  # Split Gene Name and Gene Group columns 
  library(splitstackshape)
  tmp.df <- dmps.anno.df[, .(PROBE_ID, UCSC_RefGene_Name)]
  cpg.gene.df <- cSplit(indt = tmp.df, splitCols = 'UCSC_RefGene_Name', sep = ';', type.convert = F, drop = F) %>%
    melt(id.vars = c("PROBE_ID", "UCSC_RefGene_Name"), value.name = "GeneSymbol") %>%
    dplyr::select(PROBE_ID, GeneSymbol) 
  
  # Split gene group column
  tmp.df <- dmps.anno.df[, .(PROBE_ID, UCSC_RefGene_Group)]
  cpg.gene.group.df <- cSplit(indt = tmp.df, splitCols = 'UCSC_RefGene_Group', sep = ';', type.convert = F, drop = F) %>%
    melt(id.vars = c("PROBE_ID", "UCSC_RefGene_Group"), value.name = "GeneGroup") %>%
    dplyr::select(PROBE_ID, GeneGroup) 
  
  cpg.gene.df <- cbind(cpg.gene.df, cpg.gene.group.df[, "GeneGroup"]) %>% 
    unique() # %>% na.omit()
  colnames(cpg.gene.df) <- c("PROBE_ID", "GeneSymbol", "GeneGroup")
  
  cpg.na.gene.df <- tmp.df[UCSC_RefGene_Group == "", "PROBE_ID"] %>% unique()
  
  # Get CpGs coordinates and create IRanges
  
  cpg.coord.df <- dmps.anno.df[, .(PROBE_ID, chr, pos)] %>% dplyr:: mutate(chr = sub("chr", "", chr))
  
  cpg.cord.range <-  makeGRangesFromDataFrame(cpg.coord.df, 
                                              start.field = "pos", end.field = "pos", seqnames.field = c("chr"))
  names(cpg.cord.range) <- cpg.coord.df$PROBE_ID
  
  # Get Genes' coordinates and create IRanges
  gene.list <- unique(na.omit(cpg.gene.df$GeneSymbol))
  gene.coord.df <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                         filters = 'hgnc_symbol',
                         values = gene.list,
                         mart = ensembl) %>%
    unique()
  
  gene.coord.df <- gene.coord.df[gene.coord.df$chromosome_name %in% c(1:22, "X", "Y", "MT"),]
  colnames(gene.coord.df) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")
  
  gene.cord.range <- makeGRangesFromDataFrame(gene.coord.df, 
                                              start.field = "GeneStartPos", end.field = "GeneEndPos", seqnames.field = c("GeneChr"))
  names(gene.cord.range) <- gene.coord.df$GeneSymbol
  gene.cord.range
  
  # Get Distance from CpG to all annotated genes
  
  library("gUtils")
  r1 <- cpg.cord.range
  r2 <- gene.cord.range

  dist.to.all <- gr.dist(r1, r2)
  rownames(dist.to.all) <- names(r1)
  colnames(dist.to.all) <- names(r2)
  
  dist.to.all[1:30, 1:10]
  
  dist.to.all.melt <- dist.to.all %>% 
    reshape2::melt() %>% 
    na.omit() %>% 
    unique()
  colnames(dist.to.all.melt) <- c("PROBE_ID", "GeneSymbol", "CG_GENE_DIST")
  
  dist.to.all.melt[dist.to.all.melt$PROBE_ID == "cg26056577", ]
  
  # CpG Gene Cord DF
  
  cpg.gene.coord.df <- left_join(cpg.gene.df, cpg.coord.df, by = "PROBE_ID")
  cpg.gene.coord.df <- left_join(cpg.gene.coord.df, gene.coord.df, by = "GeneSymbol")
  cpg.gene.coord.df <- cpg.gene.coord.df %>% 
    dplyr::select( PROBE_ID,  chr, pos, GeneSymbol, GeneChr, GeneStartPos, GeneEndPos, GeneGroup)
  
  cpg.gene.coord.df <- left_join(cpg.gene.coord.df, dist.to.all.melt, by = c("PROBE_ID", "GeneSymbol"))
  # cpg.gene.coord.df <- cpg.gene.coord.df[order(na.omit(cpg.gene.coord.df$CG_GENE_DIST)),] %>% setDT()
  cpg.gene.coord.df <- cpg.gene.coord.df[!is.na(cpg.gene.coord.df$CG_GENE_DIST),] %>% setDT()
  
  return(cpg.gene.coord.df)
  
}

library(dplyr)
library(reshape2)
library(IRanges)
library(ggplot2)
library(biomaRt)
ensembl <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
cpg.gene.coord.df <- GetDistances(dmps.sign.anno.df, ensembl)

cpg.gene.dist.df <- cpg.gene.coord.df[, .(PROBE_ID, GeneSymbol, CG_GENE_DIST)]
plt.df <- cpg.gene.dist.df[ , .SD[which.min(CG_GENE_DIST)], by = PROBE_ID]

# Check
plt.df[plt.df$PROBE_ID == "cg26489413", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26489413", ]

plt.df[plt.df$PROBE_ID == "cg26056577", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26056577", ]

# Join the GeneGroup
plt.df2 <- left_join(plt.df, cpg.gene.coord.df[, .(PROBE_ID, GeneSymbol, GeneGroup)])
plt.df2[plt.df2$PROBE_ID == "cg26056577", ]
plt.df2[plt.df2$PROBE_ID == "cg26489413", ]

cpg.na.gene.df[["GeneSymbol"]]   <- NA
cpg.na.gene.df[["GeneGroup"]]    <- NA
cpg.na.gene.df[["CG_GENE_DIST"]] <- NA

plt.df2 <- rbind(plt.df2, cpg.na.gene.df)

ggplot(plt.df2, aes(x = GeneGroup)) + 
  geom_bar(aes(y = (..count..)/sum(..count..), fill = GeneGroup), position = position_dodge()) + 
  geom_text(aes(label = scales::percent((..count..)/sum(..count..), accuracy = 0.1), y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -.5,  position = position_dodge(1), size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Relation to Island",
       y = "Frequency", 
       title = "The ratio of hypermethylated versus hypomethylated dDMPs in different genomic regions") +
  theme(legend.title = element_blank(), 
        legend.position = c(.9,.85),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 0.5)) 




# Look inside the obtained data
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26056577", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg15757266", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg10001905", ]

table(cpg.gene.coord.df[cpg.gene.coord.df$CG_GENE_DIST == 0, "GeneGroup"])
head(sort(table(cpg.gene.coord.df[cpg.gene.coord.df$CG_GENE_DIST == 0, "PROBE_ID"]), decreasing = T))
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26489413", ]

head(sort(table(cpg.gene.coord.df$PROBE_ID), decreasing = T))
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg27187848", ]







# Overlap

r1 <- cpg.cord.range
r2 <- gene.cord.range

overlap <- findOverlaps(r1, r2, select = "all")
dist.to.nearest.obj <- distanceToNearest(r1, r2)

hits <- names(r2)[subjectHits(dist.to.nearest.obj)]
names(hits) <- names(r1)[queryHits(dist.to.nearest.obj)]

dists <- dist.to.nearest.obj@elementMetadata@listData$distance

cpg.dist <- cbind(as.data.frame(r1[names(hits), ]), dists) %>% 
  mutate(CHR = seqnames, POS = start, DIST = dists) %>% 
  dplyr::select(CHR, POS, DIST) 

cpg.gene.df <- cbind(cpg.dist, as.data.frame(r2[hits, ]@seqnames), as.data.frame(r2[hits, ]@ranges))
colnames(cpg.gene.df) <- c("CHR", "POS", "DIST", "GENE_CHR", "GENE_START", "GENE_END", "GENE_WIDTH", "GENE_NAME")
cpg.gene.df[["PROBE_ID"]] <- rownames(cpg.gene.df)

cpg.gene.df[cpg.gene.df$PROBE_ID == "cg26056577", ]

# Distance to all

library("gUtils")
dist.to.all <- gr.dist(r1, r2)
rownames(dist.to.all) <- names(r1)
colnames(dist.to.all) <- names(r2)

dist.to.all[1:30, 1:10]

dist.to.all.melt <- dist.to.all %>% melt() %>% na.omit()
colnames(dist.to.all.melt) <- c("PROBE_ID", "GeneSymbol", "CG_GENE_DIST")

dist.to.all.melt[dist.to.all.melt$PROBE_ID == "cg26056577", ]

cpg.gene.coord.df <- left_join(cpg.gene.coord.df, dist.to.all.melt, by = c("PROBE_ID", "GeneSymbol"))
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26056577", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg15757266", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg10001905", ]
