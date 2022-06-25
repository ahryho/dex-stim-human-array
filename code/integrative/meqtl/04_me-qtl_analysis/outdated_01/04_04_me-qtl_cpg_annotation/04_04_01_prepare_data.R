# Prepare dta for the enrichment

source("~/bio/code/mpip/dex-stim-human-array/code/util.R", chdir = TRUE)

library(data.table)
library(dplyr)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)

library(annotatr)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/"

out.dir.pre <- paste0(out.dir.pre, "region_wise_independent_snps/")

# 1. Load data
# 

meqtl.dex.fn   <- paste0(out.dir.pre, "me-qtl_cis_indp_rw_dex_fdr_005.csv")
meqtl.veh.fn   <- paste0(out.dir.pre, "me-qtl_cis_indp_rw_veh_fdr_005.csv")
meqtl.delta.fn <- paste0(out.dir.pre, "me-qtl_cis_indp_rw_delta_fdr_005.csv")

# col.names <-  c("CpG_ID", "SNP", "beta", "t-stat", "p-value", "fdr")

ind.meqtl.dex.df   <- fread(meqtl.dex.fn)
ind.meqtl.veh.df   <- fread(meqtl.veh.fn)
ind.meqtl.delta.df <- fread(meqtl.delta.fn)

ind.meqtl.df <- rbind(ind.meqtl.dex.df, ind.meqtl.veh.df, ind.meqtl.delta.df)

# 2. Save Illumina annotation

# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# anno.epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# saveRDS(anno.epic, "~/bio/code/mpip/dex-stim-human-array/data/annotation/IlluminaHumanMethylationEPICanno-ilm10b4-hg19.rds")

anno.epic <- readRDS( "~/bio/code/mpip/dex-stim-human-array/data/annotation/IlluminaHumanMethylationEPICanno-ilm10b4-hg19.rds")

# 3. Generate GRanges. Extract and save annotation for DELTA me-QTLs CpGs

meqtls.cpg.lst     <- ind.meqtl.delta.df$CpG_ID %>% unique()
meqtls.cpg.anno.df <- data.frame(anno.epic[anno.epic$Name %in% meqtls.cpg.lst, ]) %>% setDT()

fwrite(meqtls.cpg.anno.df,
       paste0(out.dir.pre, "meqtls_cpg_annotated_delta.csv"),
       quote = F, row.names = F, sep   = "\t")


input <- meqtls.cpg.anno.df[,.(Name, chr, pos, Relation_to_Island)]
meqtls.cpg.delta.coord.gr <- GenomicRanges::GRanges(seqnames = input$chr,
                                                    ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                              end = as.numeric(as.character(input$pos))),
                                                    CpG_ID = input$Name,
                                                    Relation_to_Island = input$Relation_to_Island)
names(meqtls.cpg.delta.coord.gr) <- input$Name

saveRDS(meqtls.cpg.delta.coord.gr,
        paste0(out.dir.pre, "meqtl_delta_cpgs_gr.rds"))

# 4. Extract and save annotation for VEH me-QTLs CpGs

meqtls.cpg.lst     <- ind.meqtl.veh.df$CpG_ID %>% unique()
meqtls.cpg.anno.df <- data.frame(anno.epic[anno.epic$Name %in% meqtls.cpg.lst, ]) %>% setDT()

fwrite(meqtls.cpg.anno.df,
       paste0(out.dir.pre, "meqtls_cpg_annotated_veh.csv"),
       quote = F, row.names = F, sep   = "\t")

input <- meqtls.cpg.anno.df[,.(Name, chr, pos, Relation_to_Island)]
meqtls.cpg.veh.coord.gr   <- GenomicRanges::GRanges(seqnames = input$chr,
                                                    ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                              end = as.numeric(as.character(input$pos))),
                                                    CpG_ID = input$Name,
                                                    Relation_to_Island = input$Relation_to_Island)
names(meqtls.cpg.veh.coord.gr) <- input$Name

saveRDS(meqtls.cpg.veh.coord.gr,
        paste0(out.dir.pre, "meqtl_veh_cpgs_gr.rds"))

# 5. Extract and save annotation for DEX me-QTLs CpGs

meqtls.cpg.lst     <- ind.meqtl.dex.df$CpG_ID %>% unique()
meqtls.cpg.anno.df <- data.frame(anno.epic[anno.epic$Name %in% meqtls.cpg.lst, ]) %>% setDT()

fwrite(meqtls.cpg.anno.df,
       paste0(out.dir.pre, "meqtls_cpg_annotated_dex.csv"),
       quote = F, row.names = F, sep   = "\t")

input <- meqtls.cpg.anno.df[,.(Name, chr, pos, Relation_to_Island)]
meqtls.cpg.dex.coord.gr   <- GenomicRanges::GRanges(seqnames = input$chr,
                                                    ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                              end = as.numeric(as.character(input$pos))),
                                                    CpG_ID = input$Name,
                                                    Relation_to_Island = input$Relation_to_Island)
names(meqtls.cpg.dex.coord.gr) <- input$Name

saveRDS(meqtls.cpg.dex.coord.gr,
        paste0(out.dir.pre, "meqtl_dex_cpgs_gr.rds"))

# 6. Extract and save annotation fro ALL CpGs

dex.cpg.anno.df <- fread("/Users/anastasiia_hry/bio/code/mpip/dex-stim-human-array/data/methylation/dex_cpgs_annotated.csv", 
                          select = c("Name", "chr", "pos", "Relation_to_Island"))

input <- dex.cpg.anno.df[,.(Name, chr, pos, Relation_to_Island)]
dex.cpg.all.coord.gr   <- GenomicRanges::GRanges(seqnames = input$chr,
                                                 ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                           end = as.numeric(as.character(input$pos))),
                                                 CpG_ID = input$Name,
                                                 Relation_to_Island = input$Relation_to_Island)
names(dex.cpg.all.coord.gr) <- input$Name

saveRDS(dex.cpg.all.coord.gr,
        paste0(out.dir.pre, "all_cpgs_gr.rds"))

# 7. Extract and save annotation  anno epic CpGs

input <- anno.epic[,c("Name", "chr", "pos", "Relation_to_Island")]
dex.cpg.all.coord.gr   <- GenomicRanges::GRanges(seqnames = input$chr,
                                                 ranges = IRanges::IRanges(start = as.numeric(as.character(input$pos)),
                                                                           end = as.numeric(as.character(input$pos))),
                                                 CpG_ID = input$Name,
                                                 Relation_to_Island = input$Relation_to_Island)
names(dex.cpg.all.coord.gr) <- input$Name

saveRDS(dex.cpg.all.coord.gr,
        paste0(out.dir.pre, "anno_epic_gr.rds"))

# 7. Annotation from UCSC Genome Browser

library(annotatr)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

gene.annots <- c('hg19_basicgenes', 'hg19_genes_intergenic', 'hg19_genes_intronexonboundaries')

annotations <- build_annotations(
  genome = "hg19",
  annotations = gene.annots
)

## Prepare methyl ranges

meqtls.cpg.lst     <- ind.meqtl.df$CpG_ID %>% unique()
meqtls.cpg.anno.df <- data.frame(anno.epic[anno.epic$Name %in% meqtls.cpg.lst, ]) %>% setDT()

cpg.coord.df    <- meqtls.cpg.anno.df[, .(Name, chr, pos)]# %>% dplyr:: mutate(chr = sub("chr", "", chr))
cpg.coord.range <- makeGRangesFromDataFrame(cpg.coord.df, 
                                            start.field = "pos", 
                                            end.field = "pos", 
                                            seqnames.field = c("chr"))
names(cpg.coord.range) <- cpg.coord.df$Name

## Annotate

meqtls.cpg.annotated <- annotate_regions(
  regions = cpg.coord.range,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# print(delta.meqtls.cpg.annotated)
# plot_annotation(annotated_regions = delta.meqtls.cpg.annotated)

meqtls.cpg.annotated.df <- data.frame(meqtls.cpg.annotated) %>% setDT()

## 7.2 Annotation with ChipSeeker

# ALL
# 
all.cpg.anno <- annotatePeak(unique(dex.cpg.all.coord.gr), 
                               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                               annoDb = "org.Hs.eg.db")

saveRDS(all.cpg.anno,
        paste0(out.dir.pre, "cpg_annotated_withChIPseeker_all.rds"))

# DELTA
# 
delta.cpg.anno <- annotatePeak(unique(meqtls.cpg.delta.coord.gr), 
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             annoDb = "org.Hs.eg.db")

saveRDS(delta.cpg.anno,
        paste0(out.dir.pre, "meqtls_cpg_annotated_withChIPseeker_delta.rds"))

# VEH
# 
veh.cpg.anno <- annotatePeak(unique(meqtls.cpg.veh.coord.gr), 
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             annoDb = "org.Hs.eg.db")

saveRDS(veh.cpg.anno,
        paste0(out.dir.pre, "meqtls_cpg_annotated_withChIPseeker_veh.rds"))

# DEX
# 
dex.cpg.anno <- annotatePeak(unique(meqtls.cpg.dex.coord.gr), 
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             annoDb = "org.Hs.eg.db")

saveRDS(dex.cpg.anno,
        paste0(out.dir.pre, "meqtls_cpg_annotated_withChIPseeker_dex.rds"))

## 8. Chromatin annotation

chromhmm.all.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_all_states.Rds")
colnames(values(chromhmm.all.states))[1] <- "type"

meqtls.cpg.chromhmm.annotated <- annotate_regions(
  regions = cpg.coord.range,
  annotations = chromhmm.all.states,
  ignore.strand = TRUE,
  quiet = FALSE)

saveRDS(meqtls.cpg.chromhmm.annotated,
        paste0(out.dir.pre, "meqtls_cpg_chromhmm_annotated_all_models.rds"))

# meqtls.cpg.chromhmm.annotated.df <- data.frame(meqtls.cpg.chromhmm.annotated) %>% setDT()
# meqtls.cpg.chromhmm.annotated.df[["CpG_ID"]] <- names(meqtls.cpg.chromhmm.annotated)
# 
# ### delta
# 
# delta.meqtls.cpg.chromhmm.annotated.df <- meqtls.cpg.chromhmm.annotated.df[CpG_ID %in% ind.meqtl.delta.df$CpG_ID,
#                                                                            .(CpG_ID, annot.type, annot.code)] %>% unique()
# 
# fwrite(delta.meqtls.cpg.chromhmm.annotated.df,
#        paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/",
#               "meqtls_cpg_chromhmm_annotated_delta.csv"),
#        quote = F, row.names = F, sep = "\t")
# 
# ### baseline
# 
# veh.meqtls.cpg.chromhmm.annotated.df <- meqtls.cpg.chromhmm.annotated.df[CpG_ID %in% ind.meqtl.veh.df$CpG_ID,
#                                                                          .(CpG_ID, annot.type, annot.code)] %>% unique()
# 
# fwrite(veh.meqtls.cpg.chromhmm.annotated.df,
#        paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/",
#               "meqtls_cpg_chromhmm_annotated_veh.csv"),
#        quote = F, row.names = F, sep = "\t")
# 
# ### dex
# 
# dex.meqtls.cpg.chromhmm.annotated.df <- meqtls.cpg.chromhmm.annotated.df[CpG_ID %in% ind.meqtl.dex.df$CpG_ID,
#                                                                          .(CpG_ID, annot.type, annot.code)] %>% unique()
# 
# fwrite(dex.meqtls.cpg.chromhmm.annotated.df,
#        paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/",
#               "meqtls_cpg_chromhmm_annotated_dex.csv"),
#        quote = F, row.names = F, sep = "\t")
