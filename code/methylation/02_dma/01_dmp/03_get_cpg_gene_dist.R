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
  
  gene.coord.df <- gene.coord.df[gene.coord.df$chromosome_name %in% c(1:22, "X", "Y", "MT"),]
  colnames(gene.coord.df) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")
  
  # CpG Gene Cord DF
  
  cpg.gene.coord.df <- left_join(cpg.gene.df, cpg.coord.df, by = "PROBE_ID")
  cpg.gene.coord.df <- left_join(cpg.gene.coord.df, gene.coord.df, by = "GeneSymbol")
  cpg.gene.coord.df[["CG_GENE_DIST"]] <- ifelse(cpg.gene.coord.df$GeneStartPos > cpg.gene.coord.df$pos,
                                                cpg.gene.coord.df$GeneStartPos - cpg.gene.coord.df$pos,
                                                ifelse(cpg.gene.coord.df$pos > cpg.gene.coord.df$GeneEndPos,
                                                       cpg.gene.coord.df$pos - cpg.gene.coord.df$GeneEndPos - 1,
                                                       0))
  
  cpg.gene.coord.df <- cpg.gene.coord.df %>% 
    dplyr::select( PROBE_ID,  chr, pos, GeneSymbol, GeneChr, GeneStartPos, GeneEndPos, GeneGroup, CG_GENE_DIST)
  
  cpg.gene.coord.df <- cpg.gene.coord.df[!is.na(cpg.gene.coord.df$CG_GENE_DIST),] %>% setDT()
  
  cpg.na.gene.df    <- tmp.df[UCSC_RefGene_Group == "", "PROBE_ID"] %>% unique()
  
  cpg.na.gene.df[["chr"]]   <- NA
  cpg.na.gene.df[["pos"]]   <- NA
  cpg.na.gene.df[["GeneSymbol"]]   <- NA
  cpg.na.gene.df[["GeneChr"]]   <- NA
  cpg.na.gene.df[["GeneStartPos"]]   <- NA
  cpg.na.gene.df[["GeneEndPos"]]   <- NA
  cpg.na.gene.df[["GeneGroup"]]    <- NA
  cpg.na.gene.df[["CG_GENE_DIST"]] <- NA
  
  cpg.gene.coord.df <- rbind(cpg.gene.coord.df, cpg.na.gene.df)
  
  return(cpg.gene.coord.df)
}

library(dplyr)
library(reshape2)
library(IRanges)
library(ggplot2)
library(biomaRt)

ensembl           <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
cpg.gene.coord.df <- GetDistances(dmps.sign.anno.df, ensembl)

cpg.na.gene.df <- cpg.gene.coord.df[is.na(cpg.gene.coord.df$CG_GENE_DIST),] %>% unique()
