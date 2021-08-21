
source("~/bio/code/mpip/dex-stim-human-array/code/util.R", chdir = TRUE)

pkg.list <- c("tidyverse", "dplyr", "data.table")
biocmanager.pkg.list <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# Set up params

dmps.fn      <- '~/bio/code/mpip/dex-stim-human-array/output/data/methylation/01_lmem_dnam/dnam_lmem_svs_pcs_rslt.txt'
dmps.anno.fn <- '~/bio/code/mpip/dex-stim-human-array/output/data/methylation/02_dmp/dmps_annotated.csv'
 
# Load data
dmps.df <- fread(dmps.fn, sep = "\t")

# dmps.df <- fread(dmps.anno.fn)
# dmps.anno.df <- dmps.df[,1:34]

# Prepare annotation tbl

anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmps.df$PROBE_ID, ]) #c(1:4,12:19,22:43)])

dmps.anno.df <- anno.epic.sub

# Save
write.table(dmps.anno.df, file = dmps.anno.fn, sep = "\t", row.names = F)
