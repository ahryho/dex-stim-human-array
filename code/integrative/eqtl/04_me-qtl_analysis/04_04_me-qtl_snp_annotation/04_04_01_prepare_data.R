
source("~/bio/code/mpip/dex-stim-human-array/code/util.R", chdir = TRUE)

library(data.table)
library(dplyr)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

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

# 3. Extract me-QTLs SNPs, make GRanges and save 

meqtls.snp.lst            <- ind.meqtl.delta.df$SNP %>% unique()
meqtls.snp.delta.coord.df <- snp.loc[SNP %in% meqtls.snp.lst, ]
# meqtls.snp.delta.coord.df <- left_join(ind.meqtl.delta.df, snp.loc)

input <- meqtls.snp.delta.coord.df

delta.snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                                ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                          end = as.numeric(as.character(input$pos))))
names(delta.snp.coord.range) <- meqtls.snp.delta.coord.df$SNP

# Annotate

delta.snp.anno <- annotatePeak(unique(delta.snp.coord.range), 
                               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                               annoDb = "org.Hs.eg.db")
saveRDS(delta.snp.anno,
       paste0(out.dir.pre, "meqtls_snp_annotated_withChIPseeker_delta.rds"))

# baseline

meqtls.snp.lst            <- ind.meqtl.veh.df$SNP %>% unique()
meqtls.snp.veh.coord.df   <- snp.loc[SNP %in% meqtls.snp.lst, ]
# meqtls.snp.delta.coord.df <- left_join(ind.meqtl.delta.df, snp.loc)

input <- meqtls.snp.veh.coord.df

veh.snp.coord.range <- GenomicRanges::GRanges(seqnames = paste0("chr", input$chr),
                                                ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                          end = as.numeric(as.character(input$pos))))
names(veh.snp.coord.range) <- meqtls.snp.veh.coord.df$SNP

# Annotate

veh.snp.anno <- annotatePeak(unique(veh.snp.coord.range), 
                               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                               annoDb = "org.Hs.eg.db")

saveRDS(veh.snp.anno,
        paste0(out.dir.pre, "meqtls_snp_annotated_withChIPseeker_veh.rds"))


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
