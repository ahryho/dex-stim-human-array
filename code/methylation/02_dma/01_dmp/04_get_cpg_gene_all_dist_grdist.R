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

  # Get Genes' coordinates and create IRanges
  gene.list <- unique(na.omit(cpg.gene.df$GeneSymbol))
  gene.coord.df <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                         filters = 'hgnc_symbol',
                         values = gene.list,
                         mart = ensembl) %>%
    unique()
  
  chr.list <- cpg.coord.df$chr %>% unique()
  gene.coord.df <- gene.coord.df[gene.coord.df$chromosome_name %in% chr.list,]
  colnames(gene.coord.df) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")
  
  # Ranges
  
  # Get CpGs coordinates and create IRanges
  
  cpg.coord.df <- dmps.anno.df[, .(PROBE_ID, chr, pos)] %>% dplyr:: mutate(chr = sub("chr", "", chr))
  cpg.coord.range <-  makeGRangesFromDataFrame(cpg.coord.df, 
                                              start.field = "pos", end.field = "pos", seqnames.field = c("chr"))
  names(cpg.coord.range) <- cpg.coord.df$PROBE_ID
  
  gene.coord.range <- makeGRangesFromDataFrame(gene.coord.df, 
                                              start.field = "GeneStartPos", end.field = "GeneEndPos", seqnames.field = c("GeneChr"))
  names(gene.coord.range) <- gene.coord.df$GeneSymbol
  
  
  no.cores <- detectCores() - 1
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  
  # res <- foreach(chr = 1, .combine = rbind, .packages = c('gUtils', 'dplyr')) %dopar% {
  res <- foreach(chr = chr.list, .combine = rbind, .packages = c('gUtils', 'dplyr')) %dopar% {

    r1 <- cpg.coord.range[cpg.coord.range@seqnames == chr]
    r2 <- gene.coord.range[gene.coord.range@seqnames == chr]
    
    dist.to.all <- gr.dist(r1, r2)
    rownames(dist.to.all) <- names(r1)
    colnames(dist.to.all) <- names(r2)
    
    dist.to.all.melt <- dist.to.all %>% 
      reshape2::melt() %>% 
      na.omit() %>% 
      unique()
    colnames(dist.to.all.melt) <- c("PROBE_ID", "GeneSymbol", "CG_GENE_DIST")
    
    dist.to.all.melt
  }
  
  stopImplicitCluster()
  
  dist.to.all.melt <- as.data.frame(res)

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

require(foreign)

library(parallel)
library(foreach)
library(doParallel)

library(gUtils)

ensembl           <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
cpg.gene.coord.df <- GetDistances(dmps.sign.anno.df, ensembl)

cpg.na.gene.df   <- cpg.gene.coord.df[is.na(CG_GENE_DIST), .(PROBE_ID, GeneSymbol, GeneGroup, CG_GENE_DIST)]
cpg.gene.dist.df <- cpg.gene.coord.df[, .(PROBE_ID, GeneSymbol, CG_GENE_DIST)]
plt.df           <- cpg.gene.dist.df[ , .SD[which.min(CG_GENE_DIST)], by = PROBE_ID, ]

# Check
plt.df[plt.df$PROBE_ID == "cg26489413", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26489413", ]

plt.df[plt.df$PROBE_ID == "cg26056577", ]
cpg.gene.coord.df[cpg.gene.coord.df$PROBE_ID == "cg26056577", ]

# Join the GeneGroup
plt.df2 <- left_join(plt.df, cpg.gene.coord.df[, .(PROBE_ID, GeneSymbol, GeneGroup)])
plt.df2[plt.df2$PROBE_ID == "cg26056577", ]
plt.df2[plt.df2$PROBE_ID == "cg26489413", ]

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

