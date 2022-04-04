require(lme4)
require(foreign)
require(data.table)

library(parallel)
library(foreach)
library(doParallel)

# 3. with age, sex, BMI, and MDD status, SV{1-3}

# 1. Load data

# lmer.res.out.fn <- "lmer_all_plus_cell_counts.txt"
# input.parameters.fn <- "input_parameters.csv"

args            <- commandArgs(T)
gex.mtrx.fn     <- as.character(args[1])
pheno.fn        <- as.character(args[3])
lmer.res.out.fn <- as.character(args[4]) # gex_residuals
treatment       <- as.character(args[2]) # dex

# gex.mtrx.fn <- paste0("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/gex_mtrx_", treatment, ".csv")
# pheno.fn    <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
# 
gex.mtrx.fn     <- paste0(gex.mtrx.fn, treatment, ".csv")
# pheno.fn        <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
# lmer.res.out.fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/gex_residuals/gex_residuals"

gex.mtrx <- fread(gex.mtrx.fn) 
gex.ids  <- data.frame(gex.mtrx[, ENSG_ID])
gex.mtrx <- gex.mtrx[, -c("ENSG_ID")]

pheno    <- read.csv2(pheno.fn) #, na.strings = "#N/A") 
pheno    <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)
pheno$Sex        <- as.factor(pheno$Sex)
pheno$Age        <- as.numeric(pheno$Age)
pheno$BMI_D1     <- as.numeric(pheno$BMI_D1)
pheno$Status     <- as.factor(pheno$Status)

# SVs:
pheno$V1         <- as.numeric(as.character(pheno$V1))
pheno$V2         <- as.numeric(as.character(pheno$V2))
pheno$V3         <- as.numeric(as.character(pheno$V3))

# PCs:
pheno$PC1        <- as.numeric(pheno$PC1)
pheno$PC2        <- as.numeric(pheno$PC2)  

# Smkoing Score
pheno$DNAm_SmokingScore <- as.numeric(pheno$DNAm_SmokingScore)
 
# Take only veh and dex sampl
samples.ids <- as.character(pheno$DNA_ID)

# 2. Making sure about samples in pheno and and betas matrix in the same order
table(colnames(gex.mtrx) %in% samples.ids)
# pheno   <- pheno[match(colnames(beta.mtrx), pheno$DNAm_ID ),]
all(samples.ids == colnames(gex.mtrx))

gex.mtrx <- as.matrix(gex.mtrx)

length(pheno$Sex)
length(gex.mtrx[1,])

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

res <- foreach(gex = 1:nrow(gex.mtrx), .combine = rbind) %dopar% { 
  lm.model <- lm(gex.mtrx[gex, ] ~ pheno$Sex + pheno$Age + pheno$BMI_D1 + pheno$Status + pheno$DNAm_SmokingScore +
                                   pheno$V1 + pheno$V2 + pheno$V3 + 
                                   pheno$PC1 + pheno$PC2)
  
  residuals(lm.model)
}

stopImplicitCluster()

res <- as.data.frame(res)
res <- cbind(gex.ids, res)
colnames(res)[1] <- "ENSG_ID"

fwrite(res, 
       paste0(lmer.res.out.fn, "_", treatment, ".csv"),
       quote = F, row.names = F, sep = ";")
