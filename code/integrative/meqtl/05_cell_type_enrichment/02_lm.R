require(lme4)
require(foreign)
require(data.table)

library(parallel)
library(foreach)
library(doParallel)

# 1. Load data

args            <- commandArgs(T)
beta.mtrx.fn    <- as.character(args[1])
pheno.fn        <- as.character(args[2])
lmer.res.out.fn <- as.character(args[3]) # dnam_residuals
treatment       <- as.character(args[4]) # dex

# beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data_no_outlier/dex_methyl_beta_combat_mtrx.rds"
# pheno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/pheno_full_for_kimono.csv"

# beta.mtrx.fn <- paste0("~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_", treatment, ".csv")
# pheno.fn <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"

beta.mtrx.fn    <- "data/integrative/matrixEQTL/methyl_beta_mtrx_"
pheno.fn        <- "data/pheno/pheno_full_for_kimono.csv"
lmer.res.out.fn <- "output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals"

beta.mtrx.fn <- paste0(beta.mtrx.fn, treatment, ".csv")
beta.mtrx    <- fread(beta.mtrx.fn) 

cpg.ids   <- data.frame(beta.mtrx[, CpG_ID])
beta.mtrx <- beta.mtrx[, -c("CpG_ID")]
beta.mtrx <- as.matrix(beta.mtrx)

pheno     <- read.csv2(pheno.fn) #, na.strings = "#N/A") 
# pheno     <- pheno[!is.na(pheno$DNAm_ID), ]
pheno     <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)
# pheno$Sex        <- as.factor(pheno$Sex)
pheno$Age        <- as.numeric(pheno$Age)
pheno$BMI_D1     <- as.numeric(pheno$BMI_D1)
# pheno$Status     <- as.factor(pheno$Status)

# PCs:
pheno$PC1  <- as.numeric(pheno$PC1)
pheno$PC2  <- as.numeric(pheno$PC2)  

# Smkoing Score
pheno$DNAm_SmokingScore <- as.numeric(pheno$DNAm_SmokingScore)

bcc.df <- read.csv2( file = "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv")
# colnames(bcc.epidish.rpc.salas)[2:13] <- paste0("salas.", colnames(bcc.epidish.rpc.salas)[2:13])

# Take only veh and dex sampl
samples.ids <- as.character(pheno$DNA_ID)

# 2. Making sure about samples in the same order for all dfs
table(colnames(beta.mtrx) %in% samples.ids)
# pheno   <- pheno[match(colnames(beta.mtrx), pheno$DNAm_ID ),]
all(samples.ids == colnames(beta.mtrx))

# Put in the same orders
bcc.df <- bcc.df[match(pheno$DNAm_ID, bcc.df$DNAm_ID,), ]

cov.df <- cbind(bcc.df[, -1], pheno[, c("Sex", "Age", "BMI_D1", "Status", "DNAm_SmokingScore", "PC1", "PC2")])

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# res <- foreach(cpg = 1:3, .combine = rbind) %dopar% {
res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind) %dopar% { # .packages = 'lme4') %dopar% {
  lm.model <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df)
  mdl.coef <- summary(lm.model)$coefficients
  pvals <- sapply(colnames(cov.df), function(x){
    if (x %in% rownames(mdl.coef)){
      return(mdl.coef[x,'Pr(>|t|)'])
    }
    else {
      return(NA)
    }}
    )
}

stopImplicitCluster()

res <- as.data.frame(res)
res <- cbind(cpg.ids, res)
colnames(res)[1] <- "CpG_ID"

fwrite(res, 
       paste0(lmer.res.out.fn, "_", treatment, ".csv"),
       quote = F, row.names = F, sep = ";")