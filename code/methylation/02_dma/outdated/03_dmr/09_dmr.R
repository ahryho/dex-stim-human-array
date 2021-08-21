# DMR
pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

# Load libraries

pkg.list <- c("BiocManager", "tidyverse", "dplyr")

biocmanager.pkg.list <- c("minfi", "ENmix", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "GenomicRanges", "methyAnalysis")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# Input params

dir.prefix  <- "/Users/anastasiia_hry/bio/datasets/methylation/20_DMA/"
dmp.fn      <- paste0(dir.prefix, "02_dmp/dmps_pval_adj_bcc_pcs.txt")
dmp.anno.fn <- paste0(dir.prefix, "02_dmp/dmp_bcc_pcs_anno.csv") 
# dmp.anno.fn <- paste0(dir.prefix, "02_dmp/dmps_significant_annotated.txt")

# out.data.dir <-paste0(dir.prefix, "03_dmr/")

# Load data
dmp.df      <- read.table(dmp.fn, sep = "\t", header = T)
dmp.anno.df <- read.table(dmp.anno.fn, header = T, sep = ";")

# Transform the DMP df into bed format : chr | start | end | p| probe

dmp.bed.df <- dmp.anno.df %>% mutate( chrom = chr,#chr = as.numeric(substr(chr, 4, length(chr))),
                                      start = pos, 
                                      end = start + 1,
                                      p = Pval_adj, 
                                      probe = Name)  %>% 
                          dplyr::select(chrom, start, end, p, probe)
dmp.bed.df <- dmp.bed.df[order(dmp.bed.df$chrom, dmp.bed.df$start),]

write.table(dmp.bed.df, paste0(dir.prefix, "02_dmp/dmp_bcc_pcs.bed"), col.names = T, row.names = F, sep = ";", quote = F)

dmp.bed.df$start <- as.numeric(dmp.bed.df$start)
dmp.bed.df$end <- as.numeric(dmp.bed.df$end)
dmp.bed.df$p <- as.numeric(dmp.bed.df$p)
write.table(dmp.bed.df[dmp.bed.df$p < 0.1,], paste0(dir.prefix, "02_dmp/pval_for_py_combp.bed"), col.names = T, row.names = F, 
            sep = "\t", quote = F)


# Find DMRs

dist.cutoff.lst <- c("1000")
bin.size.lst    <- c("310")
seed.lst        <- c("0.01")

setwd(paste0(dir.prefix, "03_dmr/"))

dmp.sign.bed.df <- read.csv2( paste0(dir.prefix, "02_dmp/dmp_sign_bcc_pcs.bed"), header = T)
dmp.bed.df <- read.csv2( paste0(dir.prefix, "02_dmp/dmp_bcc_pcs.bed"), header = T)
dmp.bed.df$p <- as.numeric(dmp.bed.df$p)

get.acf(dmp.bed.df[dmp.bed.df$p < 1,], 1000, 310)

combp(dmp.bed.df[dmp.bed.df$p < 1,], 500, 50, 0.01, region_plot = F, mht_plot = F)
dmr.df <- read.csv("resu_combp.csv")
dmr.df <- dmr.df[order(dmr.df$chr),]
dim(dmr.df[dmr.df$nprobe >= 2, ])
dmr.df[dmr.df$nprobe >= 2, ]

combp(dmp.bed.df[dmp.bed.df$p <= 0.5 & dmp.bed.df$chr == 1,], 1000, 500, 0.4, region_plot = F, mht_plot = F)
dmr.df <- read.csv("resu_combp.csv")
dmr.df <- dmr.df[order(dmr.df$chr),]
dim(dmr.df[dmr.df$nprobe >= 2, ])
dmr.df[dmr.df$nprobe >= 2, ]

ipdmr(dmp.bed.df[dmp.bed.df$p <= 0.1,], dist.cutoff = 1000, bin.size = 310, seed = 0.01, 
      region_plot = F, mht_plot = F)
dmr.df <- read.csv("resu_ipdmr.csv")
dmr.df <- dmr.df[order(dmr.df$chr),]
dim(dmr.df[dmr.df$nprobe >= 2, ])
dmr.df[dmr.df$nprobe >= 5, ]


# Annotate DMRs

# anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmp.sign.ids, c(1:4,12:19,22:43)])

dmr.irange <- makeGRangesFromDataFrame(dmr.df, 
                                      start.field = "start", 
                                      end.field = "end", 
                                      seqnames.field = c("chr"),
                                      keep.extra.columns = T)
names(dmr.irange) <- 1:nrow(dmr.df)
dmp.bed.irange <- makeGRangesFromDataFrame(dmp.bed.df, 
                                           start.field = "start", 
                                           end.field = "end", 
                                           seqnames.field = c("chr"),
                                           keep.extra.columns = T)

dmr.anno <- annotateDMRInfo(dmr.irange, "TxDb.Hsapiens.UCSC.hg19.knownGene", promoterRange = 2000, 
                            as.GRanges = TRUE, CpGInfo = dmp.bed.irange)

dmr.anno

overlap <- findOverlaps(dmp.bed.irange, dmr.irange)
hits <- names(dmr.irange)[subjectHits(overlap)]
names(hits) <- names(dmp.bed.irange)[queryHits(overlap)]
hits[hits == "1"]

cpgs.lst <- names(hits[hits == "1"])
# cg21350406 cg18107759 cg16275316 cg16322681 cg00430562 

hits[hits == "2"]
# cg03415603 cg03104421 

dmp.bed.df[dmp.bed.df$probe %in% cpgs.lst, ]
head(dmp.bed.df[order(dmp.bed.df$p),], 20)
