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

# 
#   # Get Genes' coordinates and create IRanges
#   gene.list <- unique(na.omit(cpg.gene.df$GeneSymbol))
#   gene.coord.df <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
#                          filters = 'hgnc_symbol',
#                          values = gene.list,
#                          mart = ensembl) %>%
#     unique()
#   
#   chr.list <- dmps.anno.df$chr %>% unique() %>% sub("chr", "", .)
#   gene.coord.df <- gene.coord.df[gene.coord.df$chromosome_name %in% chr.list,]
#   colnames(gene.coord.df) <- c("GeneSymbol", "GeneChr", "GeneStartPos", "GeneEndPos")
  
  # CpG Gene Cord DF
  
  cpg.gene.coord.df <- left_join(cpg.gene.df, cpg.coord.df, by = "PROBE_ID")
 
    cpg.gene.coord.df <- left_join(cpg.gene.coord.df, gene.coord.df, by = "GeneSymbol")
    cpg.gene.coord.df[["CG_GENE_DIST"]] <- ifelse(cpg.gene.coord.df$GeneStartPos > cpg.gene.coord.df$pos,
                                                  cpg.gene.coord.df$GeneStartPos - cpg.gene.coord.df$pos,
                                                  ifelse(cpg.gene.coord.df$pos > cpg.gene.coord.df$GeneEndPos,
                                                       cpg.gene.coord.df$pos - cpg.gene.coord.df$GeneEndPos - 1,
                                                       0))
  
  #   cpg.gene.coord.df[["CG_GENE_DIST"]] <- lapply(cpg.gene.coord.df, function(x)
  #   #ifelse((x["pos"] > x["GeneStartPos"] ) & (x["pos"] < x["GeneEndPos"]),
  #   #       0,
  #          min(abs(x["pos"] - x["GeneStartPos"], abs(x["pos"] - x["GeneEndPos"])))
  # )
    
    
    # ifelse((cpg.gene.coord.df$pos > cpg.gene.coord.df$GeneStartPos ) &
    #                                               (cpg.gene.coord.df$pos < cpg.gene.coord.df$GeneEndPos),
    #                                             0,
    #                                             min(abs(cpg.gene.coord.df$pos - cpg.gene.coord.df$GeneStartPos),
    #                                                 abs(cpg.gene.coord.df$pos - cpg.gene.coord.df$GeneEndPos)))

    cpg.gene.coord.df <- cpg.gene.coord.df %>% 
      dplyr::select( PROBE_ID,  chr, pos, GeneSymbol, GeneChr, GeneStartPos, GeneEndPos, GeneGroup, CG_GENE_DIST)
  
    cpg.gene.coord.df <- cpg.gene.coord.df[!is.na(cpg.gene.coord.df$CG_GENE_DIST),] # %>% setDT()
  
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
library(ggplot2)
library(biomaRt)

ensembl           <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl')
cpg.gene.coord.df <- GetDistances(dmps.anno.df, ensembl)

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



########

chr.list <- dmps.anno.df$chr %>% unique()# %>% sub("chr", "", .)

require(foreign)
library(parallel)
library(foreach)
library(doParallel)

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

distances <- foreach(chrom = chr.list[1], .combine = rbind, .packages = c('dplyr', 'reshape2', 'biomaRt')) %dopar% {
# distances <- foreach(chrom = chr.list, .combine = rbind, .packages = c('dplyr', 'reshape2', 'biomaRt')) %dopar% {
  df <- dmps.sign.anno.df[dmps.sign.anno.df$chr == chrom,]
  GetDistances(df, ensembl)
}

stopImplicitCluster()

cpg.gene.coord.df <- distances %>% setDT()
