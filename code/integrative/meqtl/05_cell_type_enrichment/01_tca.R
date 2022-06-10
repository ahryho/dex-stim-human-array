# Cell-type enrichment
# 
# Load packages
#

library(data.table)
library(dplyr) 

library(parallel)

library(TCA)

# setwd("~/bio/code/mpip/dex-stim-human-array/")
setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/")

# Set up parameters and load data
# ####

args            <- commandArgs(T)
treatment       <- as.character(args[1])
chr.i           <- as.character(args[2])
tca.reg         <- as.character(args[3])

# Beta mtrx
# ####
methyl.mtrx.fn  <- paste0("data/methylation/beta_mtrx_by_chrom/methyl_beta_mtrx_", treatment, "_chr_", chr.i, ".csv")
methyl.mtrx.veh <- fread(methyl.mtrx.fn) 

# Pheno 
# ####

pheno.full            <- fread("data/pheno/pheno_full_for_kimono.csv", dec = c(",")) 
bcc.epidish.rpc.salas <- read.csv2( file = "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv")
colnames(bcc.epidish.rpc.salas)[2:13] <- paste0("salas.", colnames(bcc.epidish.rpc.salas)[2:13])

pheno                 <- pheno.full[Include == 1][Dex == 0]
pheno                 <- left_join(pheno, bcc.epidish.rpc.salas)

cols   <- colnames(bcc.epidish.rpc.salas)[2:13]
bcc.df <- data.frame(pheno[, ..cols], row.names = pheno$DNA_ID) # W

cov.lst <- c(
  "Sex", "Age", "BMI_D1", "DNAm_SmokingScore",
  "PC1", "PC2")

if (tca.reg == F) cov.lst <- c(cov.lst, "Status")

cov.df <- data.frame(pheno[, ..cov.lst], row.names = pheno$DNA_ID)

methyl.mtrx <- data.frame(methyl.mtrx.veh, row.names = methyl.mtrx.veh$CpG_ID) %>% select(-CpG_ID)

# Run TCA
# ####

tca.mdl <- tca(
  X = methyl.mtrx,
  W = bcc.df,
  C1 = cov.df,
  # refit_W = T,
  parallel = T,
  num_cores = round(detectCores() / 2, 0),
  max_iters = 3,
  log_file = paste0("output/data/integrative/cell_type_enrichment/tca_logs/", treatment, "_chr_", chr.i, ".txt"))

if (tca.reg == F){
  tcareg.mdl <- NULL} else{
    tcareg.mdl <- tcareg(
      X = methyl.mtrx,
      tca.mdl = tca.mdl,
      y =  data.frame(pheno[, Status], row.names = pheno$DNA_ID),
      C3 = cov.df,
      test = "marginal_conditional",
      fast_mode = T,
      parallel = T,
      num_cores = round(detectCores() / 2, 0)
    )
    
  }

result <- list(tca_mdl = tca.mdl, tcareg_mdl = tcareg.mdl)
out.fn <- paste0("output/data/integrative/cell_type_enrichment/", treatment, "/dnam_cell_type_enrichment_", treatment, "_chr_", chr.i, ".RDS") # "_chr_22.RDS")
saveRDS(result, file = out.fn)
