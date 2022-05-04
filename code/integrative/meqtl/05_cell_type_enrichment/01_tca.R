# Cell-type enrichment
# 
# Load packages
#

library("TCA")
library("ggplot2")
library("ggpubr")
library("pracma")
library("matrixStats")

 
# Load data
# 
setwd("~/bio/code/mpip/dex-stim-human-array/")

# Beta mtrx
# 
methyl.mtrx.veh <- fread("data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv") 

# Pheno 
# 
# 
pheno.full            <- fread("data/pheno/pheno_full_for_kimono.csv") 
bcc.epidish.rpc.salas <- read.csv2( file = "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv")
pheno                 <- pheno.full[Include == 1][Dex == 0]
pheno                 <- left_join(pheno.veh, bcc.epidish.rpc.salas)

cols   <- c("DNA_ID", colnames(bcc.epidish.rpc.salas))
bcc.df <- pheno[, ..cols] # W

# Group CpGs
#
out.dir.pre <- "output/data/integrative/matrixEQTL/meqtls/primary_meqtl_results/"

meqtl.parallel.fc.grp.veh.dex.gr <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_veh_dex_gr.rds"))
meqtl.parallel.fc.grp.delta.gr   <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_delta_gr.rds"))
meqtl.opposite.fc.gr             <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_opposite_fc_grp_gr.rds"))

# Prepare data
# 
parallel.fc.grp.veh.dex.cpg <- list(CpG_ID = names(meqtl.parallel.fc.grp.veh.dex.gr@anno))
parallel.fc.grp.delta.cpg   <- list(CpG_ID = names(meqtl.parallel.fc.grp.delta.gr@anno))
opposite.fc.grp.cpg         <- list(CpG_ID = names(meqtl.opposite.fc.gr@anno))

# Oposite 
# 
methyl.mtrx <- methyl.mtrx.veh[CpG_ID %in% opposite.fc.grp.cpg$CpG_ID, ]


data_path <- "./" # where files will be stored
# Download data files
download.file("github.com/cozygene/TCA/blob/master/vignettes/hannum.chr22.RData?raw=true",
                paste(data_path,"hannum.chr22.RData",sep=""))
download.file("github.com/cozygene/TCA/blob/master/vignettes/liu.cd4.chr22.RData?raw=true",
              paste(data_path,"liu.cd4.chr22.RData",sep="")) 

load("hannum.chr22.RData")
load("liu.cd4.chr22.RData")

dim(hannum$X)
hannum$X[1:10, 1:10]
hannum$cov[1:10, 1:10]
hannum$W

tca.mdl.hannum <- tca(X = hannum$X,
                      W = hannum$W,
                      C1 = hannum$cov[,c("gender","age")])
