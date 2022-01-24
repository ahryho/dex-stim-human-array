
source("~/bio/code/mpip/dex-stim-human-array/code/util.R", chdir = TRUE)

library(data.table)
library(dplyr)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)

library(annotatr)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

system(paste0("ls -lha ", out.dir.pre))

# 1. Load data

meqtl.dex.fn   <- paste0(out.dir.pre, "me-qtl_cis_indp_result_dex_fdr_005.csv")
meqtl.veh.fn   <- paste0(out.dir.pre, "me-qtl_cis_indp_result_veh_fdr_005.csv")
meqtl.delta.fn <- paste0(out.dir.pre, "me-qtl_cis_indp_result_delta_fdr_005.csv")

col.names <-  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")

ind.meqtl.dex.df   <- fread(meqtl.dex.fn, col.names = col.names)
ind.meqtl.veh.df   <- fread(meqtl.veh.fn, col.names = col.names)
ind.meqtl.delta.df <- fread(meqtl.delta.fn, col.names = col.names)

ind.meqtl.df <- rbind(ind.meqtl.dex.df, ind.meqtl.veh.df, ind.meqtl.delta.df)

# 2. Load SNP locations

snp.loc <- fread("~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/snp_locations.csv")

# 3 Calculate MAF
# for the imputed snps df from Janine

system("cd ~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_janine/qc")
system("plink2 --bfile Dex_genoData_SNPs --keep ../../../dna_ids_to_keep_for_analysis.txt --freq --out dex_geno_maf_for_filtered_samples")

# 3. Extract me-QTLs SNPs, make GRanges and save 

meqtls.snp.lst            <- ind.meqtl.delta.df$SNP %>% unique()
meqtls.snp.delta.coord.df <- snp.loc[SNP %in% meqtls.snp.lst, ]
# meqtls.snp.delta.coord.df <- left_join(ind.meqtl.delta.df, snp.loc)

input <- meqtls.snp.delta.coord.df

delta.snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                                ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                          end = as.numeric(as.character(input$pos))))
names(delta.snp.coord.range) <- meqtls.snp.delta.coord.df$SNP

makeGRangesFromDataFrame(input,start.field ="pos",end.field="pos",seqnames="chr",ignore.strand=T,keep.extra.columns=F)

# Annotate

# delta 

delta.snp.anno <- annotatePeak(unique(delta.snp.coord.range), 
                               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                               annoDb = "org.Hs.eg.db")
saveRDS(delta.snp.anno,
       paste0(out.dir.pre, "meqtls_snp_annotated_withChIPseeker_delta.rds"))

# baseline

meqtls.snp.lst            <- ind.meqtl.veh.df$SNP %>% unique()
meqtls.snp.veh.coord.df   <- snp.loc[SNP %in% meqtls.snp.lst, ]

input <- meqtls.snp.veh.coord.df

veh.snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                                ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                          end = as.numeric(as.character(input$pos))))
names(veh.snp.coord.range) <- meqtls.snp.veh.coord.df$SNP

# dex

meqtls.snp.lst            <- ind.meqtl.dex.df$SNP %>% unique()
meqtls.snp.dex.coord.df   <- snp.loc[SNP %in% meqtls.snp.lst, ]
# meqtls.snp.delta.coord.df <- left_join(ind.meqtl.delta.df, snp.loc)

input <- meqtls.snp.dex.coord.df

dex.snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                              ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                        end = as.numeric(as.character(input$pos))))
names(dex.snp.coord.range) <- meqtls.snp.dex.coord.df$SNP


# Annotate

dex.snp.anno <- annotatePeak(unique(dex.snp.coord.range), 
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             annoDb = "org.Hs.eg.db")

saveRDS(dex.snp.anno,
        paste0(out.dir.pre, "meqtls_snp_annotated_withChIPseeker_dex.rds"))


# bkgr
input <- snp.loc
snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                              ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                        end = as.numeric(as.character(input$pos))))
names(veh.snp.coord.range) <- meqtls.snp.veh.coord.df$SNP

snp.anno <- annotatePeak(unique(snp.coord.range), 
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             annoDb = "org.Hs.eg.db")

saveRDS(snp.anno,
        paste0(out.dir.pre, "snp_annotated_withChIPseeker_all.rds"))

## 7. Chromatin annotation

chromhmm.all.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_all_states.Rds")
colnames(values(chromhmm.all.states))[1] <- "type"
chromhmm.all.states <- subset(chromhmm.all.states, !(seqnames %in% c("chrM", "chrX", "chrY")))

# Subset Blood and T-cells

chromhmm.epigenomes <- read.csv("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/epigenomes.tsv", sep = "\t")
chromhmm.blood.eids <- chromhmm.epigenomes[chromhmm.epigenomes$GROUP == "Blood & T-cell", "EID"]

chromhmm.blood.states <- chromhmm.all.states[(elementMetadata(chromhmm.all.states)[, "code"] %in% chromhmm.blood.eids)]

saveRDS(chromhmm.blood.states, "~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_blood_states.Rds")

chromhmm.blood.states <- subset(chromhmm.blood.states, !(seqnames %in% c("chrM", "chrX", "chrY")))
#

AnnotateChromHMM <- function(snp.coord.range, chromhmm.all.states, out){
  x <-distanceToNearest(own, chromhmm.blood.states)
  sum(overlapsAny(delta.snp.coord.range, own))
  sum(distanceToNearest(delta.snp.coord.range, chromhmm.blood.states)@elementMetadata@listData$distance == 0)
  
  table(x@elementMetadata@listData$distance == 0)
  meqtls.snp.chromhmm.annotated <- annotate_regions(
    regions = snp.coord.range,
    annotations = chromhmm.all.states,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  meqtls.snp.chromhmm.annotated.df          <- data.frame(meqtls.snp.chromhmm.annotated) %>% setDT()
  meqtls.snp.chromhmm.annotated.df[["SNP"]] <- names(meqtls.snp.chromhmm.annotated)
  meqtls.snp.chromhmm.annotated.df          <- meqtls.snp.chromhmm.annotated.df[, .(SNP, annot.type, annot.code)] %>% unique()
  
  fwrite(meqtls.snp.chromhmm.annotated.df,
         out,
         quote = F, row.names = F, sep = "\t")
  
  return(meqtls.snp.chromhmm.annotated.df)
}

# dex

out.fn <- paste0(out.dir.pre, "meqtls_snp_chromhmm_annotated_dex.csv")
dex.meqtls.snp.chromhmm.annotated.df <- AnnotateChromHMM(delta.snp.coord.range, chromhmm.all.states, out.fn)

# delta

out.fn <- paste0(out.dir.pre, "meqtls_snp_chromhmm_annotated_delta.csv")
delta.meqtls.snp.chromhmm.annotated.df <- AnnotateChromHMM(delta.snp.coord.range, chromhmm.all.states, out.fn)

# baseline

out.fn <- paste0(out.dir.pre, "meqtls_snp_chromhmm_annotated_veh.csv")
veh.meqtls.snp.chromhmm.annotated.df <- AnnotateChromHMM(veh.snp.coord.range, chromhmm.all.states, out.fn)


# Calculate odds ratio

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
    sample_non_overlap <- length(own) - sample_overlap
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

n.delta.snps <- delta.meqtls.snp.chromhmm.annotated.df$SNP %>% unique() %>% length()
n.veh.snps   <- veh.meqtls.snp.chromhmm.annotated.df$SNP %>% unique() %>% length()

own.chromhmm <- delta.meqtls.snp.chromhmm.annotated.df[annot.type == "1_TssA", ] 
overlap      <- own.chromhmm$SNP %>% unique() %>% length()
non_overlap  <- n.delta.snps - own.chromhmm$SNP %>% unique() %>% length()

background.chromhmm   <- veh.meqtls.snp.chromhmm.annotated.df[annot.type == "1_TssA", ] 

background_bins <- lapply(chromhmm.epigenomes$EID, function(x) background.chromhmm[background.chromhmm$annot.code == x])  # get all background SNPs in MAF bin 1 to 11
own_bin_lengths <- sapply(chromhmm.epigenomes$EID, function(x) nrow(own.chromhmm[own.chromhmm$annot.code == x]))  # get length of all 11 MAF bins

nperm <- 10
resampling <- lapply(1:nperm, function(x){
  sample_overlap <- sum(sapply(1:127, function(y)
    # check that the background maf bin is not empty (happens sometimes, if the overlap is empty)
    ifelse(nrow(background_bins[[y]]) != 0, 
           # take n random samples of the background (n = own bin length except background bin length is smaller, than take all background genes of this MAF bin) 
           sample_n(background_bins[[y]], 
                    ifelse(nrow(background_bins[[y]]) >= own_bin_lengths[y], own_bin_lengths[y], nrow(background_bins[[y]]))) %>% nrow, 
           0)
  )
  )
  sample_non_overlap <- n.veh.snps - sample_overlap
  conf_mtrx <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
  fisher.test.rslt <- fisher.test(conf_mtrx)
  c(p_value = fisher.test.rslt$p.value, fisher.test.rslt$estimate)
}
)

p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] <= p_value) + 1) / (nperm + 1)

p_value <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
or <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
c(p_value, or)

# public       <- chromhmm.all.states[(elementMetadata(chromhmm.all.states)[, "code"] %in% "E033")]
