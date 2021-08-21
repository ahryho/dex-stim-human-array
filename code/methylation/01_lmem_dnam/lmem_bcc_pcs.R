require(lme4)
require(foreign)

library(parallel)
library(foreach)
library(doParallel)

# 3. with age, sex, BMI, and MDD status, BCCs , PC1 and PC2 as covariates

# 1. Load data

# lmer.res.out.fn <- "lmer_all_plus_cell_counts.txt"
# input.parameters.fn <- "input_parameters.csv"

args            <- commandArgs(T)
beta.mtrx.fn    <- as.character(args[1])
pheno.fn        <- as.character(args[2])
lmer.res.out.fn <- as.character(args[3])

# beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# pheno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/pheno_full_for_kimono.csv"

# beta.mtrx.fn <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# pheno.fn <- "~/bio/datasets/pheno/pheno_full_for_kimono.csv"

# fwrite(pheno, 
#        pheno.fn,
#        row.names = F,
#        quote = F,
#        sep = ";",
#        dec = ".")

beta.mtrx <- readRDS(beta.mtrx.fn) 
pheno     <- read.csv2(pheno.fn)#, na.strings = "#N/A") 
#pheno     <- pheno[!is.na(pheno$DNAm_ID), ]
pheno     <- pheno[pheno$Include == 1, ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)
pheno$Sex        <- as.factor(pheno$Sex)
# pheno$Group      <- as.factor(pheno$Group)
pheno$Age        <- as.numeric(pheno$Age)
pheno$BMI        <- as.numeric(pheno$BMI)
pheno$Status     <- as.factor(pheno$Status)

# Cell types:
pheno$CD8T  <- as.numeric(pheno$CD8T)
pheno$CD4T  <- as.numeric(pheno$CD4T)
pheno$NK    <- as.numeric(pheno$NK)
pheno$Bcell <- as.numeric(pheno$Bcell)
pheno$Mono  <- as.numeric(pheno$Mono)
pheno$Gran  <- as.numeric(pheno$Gran)

# PCs:
pheno$PC1  <- as.numeric(pheno$PC1)
pheno$PC2  <- as.numeric(pheno$PC2)  

# Take only veh and dex sampl
samples.veh.ids <- pheno$DNAm_ID[pheno$Group == "veh"]
samples.dex.ids <- pheno$DNAm_ID[pheno$Group == "dex"]

# 2. Making sure about samples in pheno and and betas matrix in the same order

table(colnames(beta.mtrx) %in% pheno$DNAm_ID)

all(pheno$DNAm_ID == colnames(beta.mtrx))

pheno   <- pheno[match(colnames(beta.mtrx), pheno$DNAm_ID ),]
all(pheno$DNAm_ID == colnames(beta.mtrx))

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# res <- foreach(cpg =  1:10, .combine = rbind, .packages = 'lme4') %dopar% {
res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind, .packages = 'lme4') %dopar% {
  
  # for (cpg in 1:nrow(beta.mtrx)){
  sample <- pheno$Sample_ID #pheno$Individual
  
  lmer.null  <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Mono + pheno$Gran +
                       pheno$PC1 + pheno$PC2 + 
                       (1|sample), REML = F)
  
  lmer.model <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Mono + pheno$Gran + 
                       pheno$PC1 + pheno$PC2 + 
                       pheno$Group + (1|sample), REML = F)
  
  res.anova  <- anova(lmer.null, lmer.model)
  
  # Statistics
  mean.veh <- round(mean(beta.mtrx[cpg, samples.veh.ids]), 4)
  mean.dex <- round(mean(beta.mtrx[cpg, samples.dex.ids]), 4)
  sd.veh   <- round(sd(beta.mtrx[cpg, samples.veh.ids]), 4)
  sd.dex   <- round(sd(beta.mtrx[cpg, samples.dex.ids]), 4)
  var.veh  <- round(var(beta.mtrx[cpg, samples.veh.ids]), 4)
  var.dex  <- round(var(beta.mtrx[cpg, samples.dex.ids]), 4)
  
  fc       <- round(mean.dex - mean.veh, 4)
  var.all  <- round(var(beta.mtrx[cpg, ]), 4)
  mad      <- round(mad(beta.mtrx[cpg, ]), 4)
  
  cbind(as.character(rownames(beta.mtrx)[cpg]), 
        signif(res.anova$"Pr(>Chisq)"[2], 4),
        round(res.anova$Chisq[2], 4),
        mean.veh, sd.veh, var.veh,
        mean.dex, sd.dex, var.dex,
        fc, var.all, mad)
}

stopImplicitCluster()

fdr <- p.adjust(res[, 2], method = "fdr")
res <- cbind(res, fdr)

res.colnames <- cbind("PROBE_ID", "PVAL", "CHI_SQ",
                      "MEAN_VEH", "SD_VEH", "VAR_VEH",
                      "MEAN_DEX", "SD_DEX", "VAR_DEX",
                      "FC", "VAR_ALL", "MAD", "FDR")
res          <- rbind(res.colnames, res)

write.table(res,
            file = lmer.res.out.fn, row.names = F, quote = F, sep = "\t", col.names = F, append = F)