library(dplyr)
library(splitstackshape)
library(reshape2)
library(IRanges)
library(ggplot2)
library(biomaRt)

require(foreign)

library(parallel)
library(foreach)
library(doParallel)

library(gUtils)

# Load data

dmps.anno.fn <- '~/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dex_cpgs_annotated.csv'
# dmps.anno.fn <- '/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dex_cpgs_annotated.csv'
dmps.anno.df <- fread(dmps.anno.fn, sep = "\t")
colnames(dmps.anno.df)[4] <- "PROBE_ID"

gene.expression.fn <- '~/bio/code/mpip/dex-stim-human-array/data/gene_expression/mapping_ilmn_ensg_gene.csv'
# gene.expression.fn <- '/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/gene_expression/mapping_ilmn_ensg_gene.csv'
gex.df <- fread(gene.expression.fn)

# Prepare methyl ranges

cpg.coord.df <- dmps.anno.df[, .(PROBE_ID, chr, pos)] %>% dplyr:: mutate(chr = sub("chr", "", chr))
cpg.coord.range <-  makeGRangesFromDataFrame(cpg.coord.df, 
                                             start.field = "pos", end.field = "pos", seqnames.field = c("chr"))
names(cpg.coord.range) <- cpg.coord.df$PROBE_ID

chr.list <- cpg.coord.df$chr %>% unique()

# Prepare gex ranges

ensembl  <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
gene.list <- unique(na.omit(gex.df$Illumina_ID))
gene.coord.df <- getBM(attributes = c('illumina_humanht_12_v3', 'chromosome_name', 'start_position', 'end_position'),
                       filters = 'illumina_humanht_12_v3',
                       values = gene.list,
                       mart = ensembl) %>%
  unique()

gene.coord.df <- gene.coord.df[gene.coord.df$chromosome_name %in% chr.list,]
colnames(gene.coord.df) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")

gene.coord.range <- makeGRangesFromDataFrame(gene.coord.df, 
                                             start.field = "GeneStartPos", end.field = "GeneEndPos", seqnames.field = c("GeneChr"))
names(gene.coord.range) <- gene.coord.df$GeneSymbol

# Find distances

# no.cores <- detectCores() - 1
# cl <- makeCluster(no.cores - 10)
# registerDoParallel(no.cores)
registerDoParallel(cores = 50)

distances <- foreach(chr = chr.list, .combine = rbind, .packages = c('gUtils', 'dplyr')) %dopar% {
  
  r1 <- cpg.coord.range[cpg.coord.range@seqnames == chr] # GRanges object with 740357 ranges and 0 metadata columns
  r2 <- gene.coord.range[gene.coord.range@seqnames == chr] # GRanges object with 12841 ranges and 0 metadata columns
  
  dist.to.all <- gr.dist(r1, r2)
  rownames(dist.to.all) <- names(r1)
  colnames(dist.to.all) <- names(r2)
  
  dist.to.all.melt <- dist.to.all %>% 
    reshape2::melt() %>% 
    na.omit() %>% 
    unique()
  colnames(dist.to.all.melt) <- c("PROBE_ID", "ILMN_ID", "CG_GENE_DIST")
  
  return(dist.to.all.melt)
}

stopImplicitCluster()

colnames(distances) <- c("PROBE_ID", "ILMN_ID", "CG_GENE_DIST")

cpg.gene.coord.df <- distances %>% setDT()

cpg.gene.coord.df <- left_join(cpg.gene.coord.df, cpg.coord.df, by = "PROBE_ID")
cpg.gene.coord.df <- left_join(cpg.gene.coord.df, unique(gene.coord.df), by = c("ILMN_ID" = "GeneSymbol"))
cpg.gene.coord.df <- cpg.gene.coord.df %>% dplyr::select(PROBE_ID, CG_CHR = chr, CG_POS = pos,
                                                  ILMN_ID, GENE_CHR = GeneChr, GENE_START_POS = GeneStartPos, GENE_END_POS = GeneEndPos,
                                                  CG_GENE_DIST)

registerDoParallel(cores = 50)

foreach(chr = chr.list, .combine = rbind, .packages = c('data.table')) %dopar% {
  fwrite(cpg.gene.coord.df[CG_CHR == chr],
         paste0('/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/dex_cpgs_ilmn_genes_distances/dex_cpgs_ilmn_genes_distances_chr_',
                chr, '.csv'),
         sep = ";", quote = F, row.names = F)
  return(chr)
}

stopImplicitCluster()

fwrite(cpg.gene.coord.df,
       '/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/dex_cpgs_ilmn_genes_distances_pe5.csv',
       sep = ";", quote = F, row.names = F)

# Save distances which are less then 1mbp
cpg.gene.coord.dist.1mbp.df <- cpg.gene.coord.df[CG_GENE_DIST <= 1000000]
fwrite(cpg.gene.coord.dist.1mbp.df,
       '/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/dex_cpgs_ilmn_genes_distances_upto_1mbp.csv',
       sep = ";", quote = F, row.names = F)

# Save distances which are less then 10 mbp
cpg.gene.coord.dist.10mbp.df <- cpg.gene.coord.df[CG_GENE_DIST <= 10000000]
fwrite(cpg.gene.coord.dist.10mbp.df,
       '/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/dex_cpgs_ilmn_genes_distances_upto_10mbp.csv',
       sep = ";", quote = F, row.names = F)