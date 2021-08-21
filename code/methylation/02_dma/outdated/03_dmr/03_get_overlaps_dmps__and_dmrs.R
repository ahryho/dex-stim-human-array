# Find overlaps between sign dmps and dmrs

pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

pkg.list <- c("tidyverse", "dplyr", "ggplot2", "Gviz", "RColorBrewer")
biocmanager.pkg.list <- c("GenomicRanges", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

src.data.pre <- "~/bio/datasets/methylation/"

# DMPs filename
dmps.fn.sufix <- "bcc_pcs_beta_10_p_5.txt"
dmps.fn       <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_significant_with_beta_stat_", dmps.fn.sufix)

# DMRs filename

dmrs.fn.sufix <- "dex.regions-t.bed"
dmrs.fn <- paste0(src.data.pre, "20_DMA/03_dmr/", "comb-p/", dmrs.fn.sufix)
dmrs.fn.sufix <- "dex.anno.hg19.bed"
dmrs.anno.fn <- paste0(src.data.pre, "20_DMA/03_dmr/", "comb-p/", dmrs.fn.sufix)

# Load DMPs

dmps.df <- read.csv(dmps.fn, sep = "\t")

# Prepare annotation tbl

anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmps.df$Probe_Id, c(1:4)])

# Merge
dmps.anno.df <- right_join(anno.epic.sub, dmps.df, by = c("Name" = "Probe_Id"))  # 1640 x 6
# dmps.anno.df <- dmps.anno.df[order(dmps.anno.df$Pval_adj),]

# Load DMRs
dmrs.df <- read.csv(dmrs.anno.fn, sep = "\t") # dmr > 0.05

# Conver DF to IRange
dmps.irange <- makeGRangesFromDataFrame(dmps.anno.df, 
                                        start.field = "pos", 
                                        end.field = "pos", 
                                        seqnames.field = c("chr"),
                                        keep.extra.columns = T)
names(dmps.irange) <- dmps.irange$Name
dmps.irange

dmrs.irange <- makeGRangesFromDataFrame(dmrs.df, 
                                        start.field = "start", 
                                        end.field = "end", 
                                        seqnames.field = c("X.chrom"),
                                        keep.extra.columns = T)
names(dmrs.irange) <- 1:nrow(dmrs.df)
dmrs.irange

# Find overlaps
overlap <- findOverlaps(dmps.irange, dmrs.irange)
hits <- names(dmrs.irange)[subjectHits(overlap)]
names(hits) <- names(dmps.irange)[queryHits(overlap)]
hits

overlap <- findOverlaps(dmrs.irange, dmps.irange)
hits <- names(dmps.irange)[subjectHits(overlap)]
names(hits) <- names(dmrs.irange)[queryHits(overlap)]
hits

# Get DMPs that are not in DMRs
dmps.unique.df <- dmps.anno.df[!(dmps.anno.df$Name %in% hits), ]

write.table(dmps.unique.df,
            file = paste0(src.data.pre , "20_DMA/03_dmr/dmps_unique_", dmps.fn.sufix), 
            row.names = F, quote = F, sep = "\t", col.names = T)

