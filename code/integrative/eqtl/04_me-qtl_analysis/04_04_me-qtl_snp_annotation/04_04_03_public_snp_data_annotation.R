library(data.table)
library(GenomicRanges)
library(tidyr)
library(arules)
library(reshape2)
library(rtracklayer)
library(foreach)
library(doSNOW)
library(S4Vectors)

# function, to calculate the SNP annotation to a public dataset with multiple categories, but with respect to the MAF bins
# @param own          Own cis-eQTL data set as GRanges object
# @param public       A public cis-eQTL data set as GRanges object
# @param nperm        Number of permutations (usually a 1000 times)
# @param background   All background SNPs of the eQTL generation overlapping with the public dataset (changes for each public dataset)
# @return             Vector including the number of overlapping SNPs (1 = overlap), the resampling p-value (2 = p_value) and the resampling odds ratio (3 = or)
AnnotatePublicData <- function(own, public, nperm, background){
  #print(paste0("Number of overlapping background SNPs: ", length(background@ranges)))  # print number of overlapping background SNPs
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  # resample nperm times
  resampling <- lapply(1:nperm, function(x)
    sum(sapply(1:11, function(y)
      # check that the background maf bin is not empty (happens sometimes, if the overlap is empty)
      ifelse(length(background_bins[[y]]) != 0, sum(
        # take n random samples of the background (n = own bin length except background bin length is smaller, than take all background genes of this MAF bin) 
        IRanges::overlapsAny(sample(background_bins[[y]], 
                                    ifelse(length(background_bins[[y]]) >= own_bin_lengths[y], own_bin_lengths[y], length(background_bins[[y]]))), 
                             public)
      ), 0))))
  # IF 0 => actually smaller than 1e^(-5/6)
  p_value <- ifelse((sum(resampling >= overlap) / nperm) == 0, paste("<=", 1/nperm, sep = " "), sum(resampling >= overlap) / nperm)
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling p-value
}

AnnotatePublicData <- function(own, public, nperm, background){
  #print(paste0("Number of overlapping background SNPs: ", length(background@ranges)))  # print number of overlapping background SNPs
  overlap     <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  non_overlap <- length(own) - overlap
  
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  
  # resample nperm times
  resampling <- lapply(1:nperm, function(x){
    sample_overlap <- sum(sapply(1:11, function(y)
      # check that the background maf bin is not empty (happens sometimes, if the overlap is empty)
      ifelse(length(background_bins[[y]]) != 0, 
             # take n random samples of the background (n = own bin length except background bin length is smaller, than take all background genes of this MAF bin) 
             sum(IRanges::overlapsAny(sample(background_bins[[y]], 
                                             ifelse(length(background_bins[[y]]) >= own_bin_lengths[y], own_bin_lengths[y], length(background_bins[[y]]))), 
                                      public)), 
             0)
      )
    )
    sample_non_overlap <- length(background) - sample_overlap
    conf_mtrx <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
    fisher.test.rslt <- fisher.test(conf_mtrx)
    c(p_value = fisher.test.rslt$p.value, fisher.test.rslt$estimate)
  }
  )
  # IF 0 => actually smaller than 1e^(-5/6)
  p_value <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  return(c(p_value, or)) # return the resampling p-value
}


out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

background.gr      <- readRDS(paste0(out.dir.pre, "meqtl_veh_snps_with_maf_gr.rds"))
meqtl.delta.snp.gr <- readRDS(paste0(out.dir.pre, "meqtl_delta_snps_with_maf_gr.rds"))
gwas.gr            <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds")

# filter GWAS datasets for genome-wide (gw) hits (p-value <= 5e-8)  x

gwas <- GRangesList("CD_2019" = gwas.gr)
gwas.gw <- GRangesList("CD_2019" = gwas.gw.gr)


### Annotation to genome-wide GWAS hits ###
# First, set parameters
perm <- 1
# Annotation to GWAS summary statistics (parallelized)
cl <- makeCluster(2)
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
results2 <- foreach (i = seq_along(gwas), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                       results <- data.frame()
                       overlap <- paste0(names(gwas)[i], "_SNP_overlap")
                       pval <- paste0(names(gwas)[i], "_SNP_pVal")
                       or <- paste0(names(gwas)[i], "_SNP_or")
                       gwas_selected <- gwas[[i]]
                       background <- subsetByOverlaps(background.all.gr, gwas_selected)
                       results[c(overlap, pval, or), "Delta"] <- AnnotatePublicData(meqtl.delta.snp.gr, gwas_selected, perm, background.gr)
                       
                       # write.table(results, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas_gw)[i],"_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                       results
                     }
stopCluster(cl)
own <- meqtl.delta.snp.gr
public <- gwas_selected
nperm <- 5

results2
# Write results into a file
write.table(results2, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
rm(results, results2)