setwd("~/bio/code/mpip/dex-stim-human-array/")

# Load functions

source("code/integrative/util.R")
source("code/integrative/meqtl/05_cell_type_enrichment/04_overlap_with_provencal_pnas_data/00_custom_functions.R")

# Load libraries

library(data.table)

# Load public data

public.data.fn <- "data/public_data/Provencal2019_Hipocampal_CpGs/pnas_provencal_2019_hipocampal_cpgs_list.xlsx"
public.data    <- read.xlsx2(public.data.fn, sheetName = "S7",  colIndex = 1, header = T)

# Load  pheno data and bcc

pheno     <- LoadPheno()
bcc.df    <- LoadBCC()

# Extract covariates and combine with BBCs

cov.lst <- c("Dex", "Sex", "Age", "BMI_D1", "Status", "DNAm_SmokingScore", 
             "DNAm_SV1", "DNAm_SV2", "DNAm_SV3", "PC1", "PC2")

cov.df  <- pheno[, c("Sample_ID", "DNA_ID", "RNA_ID", "DNAm_ID", cov.lst)]

# Extract CpGs

beta.veh.sub.mtrx <- GetMethylSubset(treatment = "veh", public.data, is.out = T) # 442 x 196 
beta.dex.sub.mtrx <- GetMethylSubset(treatment = "dex", public.data, is.out = T) # 442 x 196

# Extract cell type-specific analysis results (FDR)

fdr.veh.bcc.df <- GetFDRSubset(treatment = "veh", public.data)
fdr.dex.bcc.df <- GetFDRSubset(treatment = "dex", public.data)

# Calculate EpiStress score 
# 
public.data.fn <- "data/public_data/Provencal2019_Hipocampal_CpGs/pnas_provencal_2019_hipocampal_24cpgs_with_weights.xlsx"
public.data    <- read.xlsx2(public.data.fn, sheetName = "S8", header = T)
public.data$Weights <- as.numeric(public.data$Weights)

beta.veh.sub.mtrx   <- GetMethylSubset(treatment = "veh", public.data, is.out = F) # 23  x 196 
beta.dex.sub.mtrx   <- GetMethylSubset(treatment = "dex", public.data, is.out = F) # 23 x 196
beta.delta.sub.mtrx <- GetMethylSubset(treatment = "delta", public.data, is.out = F) # 23 x 196

beta.delta.sub.mtrx <- (beta.dex.sub.mtrx[, -1] - beta.veh.sub.mtrx[, -1]) 
beta.delta.sub.mtrx <- cbind(beta.dex.sub.mtrx[, 1], beta.delta.sub.mtrx)

epistress.score.veh.df   <- GetEpiStressScore(public.data, beta.veh.sub.mtrx, "veh")
epistress.score.dex.df   <- GetEpiStressScore(public.data, beta.dex.sub.mtrx, "dex")
epistress.score.delta.df <- GetEpiStressScore(public.data, beta.delta.sub.mtrx, "delta")
