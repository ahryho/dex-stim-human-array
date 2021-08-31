
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

# # Prepare annotation tbl
# 
# anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# location.epic <- getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)   
# 
# location.epic <- location.epic[location.epic$Name %in% dmps.anno.df$PROBE_ID, ]

# Get CpGs coordinates and create IRanges

cpg.cord.tbl <- dmps.sign.anno.df[, .(PROBE_ID, chr, pos)] %>% dplyr:: mutate(chr = sub("chr", "", chr))

# Get Genes' coordinates and create IRanges

library(splitstackshape)
tmp.df <- dmps.sign.anno.df[, .(PROBE_ID, UCSC_RefGene_Name)]

cpg.gene.df <- cSplit(indt = tmp.df, splitCols = 'UCSC_RefGene_Name', sep = ';', type.convert = F, drop = F) %>%
  melt(id.vars = c("PROBE_ID", "UCSC_RefGene_Name"), value.name = "GeneSymbol") %>%
  dplyr::select(PROBE_ID, GeneSymbol) %>%
  unique()  # %>% na.omit()

gene.list <- unique(na.omit(cpg.gene.df$GeneSymbol))
library(biomaRt)
ensembl      <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
gene.cord.tbl <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                       filters = 'hgnc_symbol',
                       values = gene.list,
                       mart = ensembl) %>%
  unique()

gene.cord.tbl <- gene.cord.tbl[gene.cord.tbl$chromosome_name %in% c(1:22, "X", "Y", "MT"),]
colnames(gene.cord.tbl) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")

# CpG Gene Cord DF

cpg.gene.coord.df <- left_join(cpg.gene.df, cpg.cord.tbl, by = "PROBE_ID")
cpg.gene.coord.df <- left_join(cpg.gene.coord.df, gene.cord.tbl, by = "GeneSymbol")
cpg.gene.coord.df <- cpg.gene.coord.df %>% 
  dplyr::select( PROBE_ID,  chr, pos, GeneSymbol, GeneChr, GeneStartPos, GeneEndPos)

# Create IRanges

cpg.cord.range <-  makeGRangesFromDataFrame(cpg.cord.tbl, 
                                            start.field = "pos", end.field = "pos", seqnames.field = c("chr"))
names(cpg.cord.range) <- cpg.cord.tbl$PROBE_ID
cpg.cord.range

gene.cord.range <- makeGRangesFromDataFrame(gene.cord.tbl, 
                                            start.field = "GeneStartPos", end.field = "GeneEndPos", seqnames.field = c("GeneChr"))
names(gene.cord.range) <- gene.cord.tbl$GeneSymbol


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
