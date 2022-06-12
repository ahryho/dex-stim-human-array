require(foreign)
require(data.table)

library(parallel)
library(foreach)
library(doParallel)

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/")

# 1. Load data

args            <- commandArgs(T)
treatment       <- as.character(args[1]) # dex
# beta.mtrx.fn    <- as.character(args[2])
# pheno.fn        <- as.character(args[3])

beta.mtrx.fn    <- paste0("data/integrative/matrixEQTL/methyl_beta_mtrx_", treatment, ".csv")
pheno.fn        <- "data/pheno/pheno_full_for_kimono.csv"
bcc.fn          <- "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv" # Salas BCCs
out.fn          <- paste0("output/data/integrative/cell_type_enrichment/dnam_cell_type_enrichment_pvals_", treatment, ".csv")

# Load and prepare methyl beta mtrx

beta.mtrx       <- fread(beta.mtrx.fn) 
cpg.ids         <- data.frame(beta.mtrx[, CpG_ID])
beta.mtrx       <- beta.mtrx[, -c("CpG_ID")]
beta.mtrx       <- as.matrix(beta.mtrx)

# Load and prepare pheno data

pheno           <- read.csv2(pheno.fn) #, na.strings = "#N/A") 
dnam.veh.ids    <- pheno[pheno$Include == 1 & pheno$Group == "veh", "DNAm_ID"]
pheno           <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

pheno$Sample_ID         <- as.factor(pheno$Sample_ID)
pheno$Age               <- as.numeric(pheno$Age)
pheno$BMI_D1            <- as.numeric(pheno$BMI_D1)
pheno$PC1               <- as.numeric(pheno$PC1)
pheno$PC2               <- as.numeric(pheno$PC2)  
pheno$DNAm_SmokingScore <- as.numeric(pheno$DNAm_SmokingScore)

# BCC

bcc.df <- read.csv2(bcc.fn)
bcc.df <- bcc.df[match(dnam.veh.ids, bcc.df$DNAm_ID,), ] # put in the same order bcc df as pheno df
cov.df <- cbind(bcc.df[, -1], pheno[, c("Sex", "Age", "BMI_D1", "Status", "DNAm_SmokingScore", "PC1", "PC2")])

# 2. Making sure about samples in the same order for all dfs

samples.ids <- as.character(pheno$DNA_ID)
table(colnames(beta.mtrx) %in% samples.ids)
all(samples.ids == colnames(beta.mtrx))

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# res <- foreach(cpg = 1:3, .combine = rbind) %dopar% {
res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind) %dopar% {
  lm.model <- lm(beta.mtrx[cpg, ] ~ (. - 1), data = cov.df)
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

# 4. Save results

fwrite(res, 
       out.fn,
       quote = F, row.names = F, sep = ";")