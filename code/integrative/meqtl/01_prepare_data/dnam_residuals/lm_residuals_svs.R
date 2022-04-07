require(lme4)
require(foreign)
require(data.table)

library(parallel)
library(foreach)
library(doParallel)

# 1. Load data

# lmer.res.out.fn <- "lmer_all_plus_cell_counts.txt"
# input.parameters.fn <- "input_parameters.csv"

args            <- commandArgs(T)
beta.mtrx.fn    <- as.character(args[1])
pheno.fn        <- as.character(args[2])
lmer.res.out.fn <- as.character(args[3]) # dnam_residuals
treatment       <- as.character(args[4]) # dex

# beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data_no_outlier/dex_methyl_beta_combat_mtrx.rds"
# pheno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/pheno_full_for_kimono.csv"

beta.mtrx.fn <- paste0(beta.mtrx.fn, treatment, ".csv")
beta.mtrx <- fread(beta.mtrx.fn) 
cpg.ids   <- data.frame(beta.mtrx[, CpG_ID])
beta.mtrx <- beta.mtrx[, -c("CpG_ID")]

pheno     <- read.csv2(pheno.fn) #, na.strings = "#N/A") 
# pheno     <- pheno[!is.na(pheno$DNAm_ID), ]
pheno     <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)

# SVs:
pheno$DNAm_SV1  <- as.numeric(as.character(pheno$DNAm_SV1))
pheno$DNAm_SV2  <- as.numeric(as.character(pheno$DNAm_SV2), dec = ".")
pheno$DNAm_SV3  <- as.numeric(as.character(pheno$DNAm_SV3))

# Take only veh and dex sampl
samples.ids <- as.character(pheno$DNA_ID)

# 2. Making sure about samples in pheno and and betas matrix in the same order
table(colnames(beta.mtrx) %in% samples.ids)
# pheno   <- pheno[match(colnames(beta.mtrx), pheno$DNAm_ID ),]
all(samples.ids == colnames(beta.mtrx))

beta.mtrx <- as.matrix(beta.mtrx)

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# res <- foreach(cpg = 1:3, .combine = rbind) %dopar% {
res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind) %dopar% { # .packages = 'lme4') %dopar% {
  lm.model <- lm(beta.mtrx[cpg, ] ~ pheno$DNAm_SV1 + pheno$DNAm_SV2 + pheno$DNAm_SV3)

  residuals(lm.model)
}

stopImplicitCluster()

res <- as.data.frame(res)
res <- cbind(cpg.ids, res)
colnames(res)[1] <- "CpG_ID"

fwrite(res, 
       paste0(lmer.res.out.fn, "_", treatment, ".csv"),
       quote = F, row.names = F, sep = ";")
