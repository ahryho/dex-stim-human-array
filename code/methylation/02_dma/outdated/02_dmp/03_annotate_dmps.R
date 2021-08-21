
pre <- "~/bio/code/mpip/dex-stim-human-differential-methyl-analysis/"
setwd(pre)
source("util.R", chdir = TRUE)

pkg.list <- c("tidyverse", "dplyr")
biocmanager.pkg.list <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# Set up params

src.data.pre <- "/Users/anastasiia_hry/bio/datasets/methylation/" 
dmps.fn      <- paste0(src.data.pre, "20_DMA/01_lme_models/lme_dex_svs.txt")
dmps.anno.fn <- paste0(src.data.pre, "20_DMA/02_dmp/dmps_annotated.csv")
 
# Load data
dmps.df <- read.csv(dmps.fn, sep = "\t")

# dmps.df <- read.csv(dmps.anno.fn, sep = "\t")
# dmps.anno.df <- dmps.df[,1:34]

# Prepare annotation tbl

anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmps.df$Probe_Id, ]) #c(1:4,12:19,22:43)])

# Merge
# dmps.anno.df <- right_join(anno.epic.sub, dmps.df, by = c("Name" = "Probe_Id"))  # 740,357 x 36
# dmps.anno.df <- dmps.anno.df[order(dmps.anno.df$Pval_adj),]

dmps.anno.df <- anno.epic.sub
# Save

write.table(dmps.anno.df, file = dmps.anno.fn, sep = "\t", row.names = F)
