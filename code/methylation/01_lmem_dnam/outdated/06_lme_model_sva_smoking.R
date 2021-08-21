require(lme4)
require(foreign)

library(parallel)
library(foreach)
library(doParallel)

# Calculate lmem with age, sex, BMI, and MDD status, SV{1-3} and DNAm smoking score

# 1. Load data

# lmer.res.out.fn <- "lmer_all_plus_cell_counts.txt"
# input.parameters.fn <- "input_parameters.csv"

args            <- commandArgs(T)
beta.mtrx.fn    <- as.character(args[1])
pheno.fn        <- as.character(args[2])
lmer.res.out.fn <- as.character(args[3])

beta.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
pheno.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/pheno_full_for_kimono.csv"

# beta.mtrx.fn <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
# pheno.fn <- "~/bio/datasets/pheno/pheno_full_for_kimono.csv"

beta.mtrx <- readRDS(beta.mtrx.fn) 
pheno     <- read.csv2(pheno.fn, na.strings = "NA") 
pheno     <- pheno[!is.na(pheno$DNAm_ID), ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)
pheno$Sex        <- as.factor(pheno$Sex)
# pheno$Group      <- as.factor(pheno$Group)
pheno$Age        <- as.numeric(pheno$Age)
pheno$BMI        <- as.numeric(pheno$BMI)
pheno$Status     <- as.factor(pheno$Status)

# SVs:
pheno$DNAm_SV1  <- as.numeric(pheno$DNAm_SV1)
pheno$DNAm_SV2  <- as.numeric(pheno$DNAm_SV2)
pheno$DNAm_SV3  <- as.numeric(pheno$DNAm_SV3)

# Elliot DNAm Smoking factor
pheno$DNAm_SmokingScore <- as.numeric(pheno$DNAm_SmokingScore)


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

# res <- foreach(cpg =  1:2, .combine = rbind, .packages = 'lme4') %dopar% {
res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind, .packages = 'lme4') %dopar% {
  
  # for (cpg in 1:nrow(beta.mtrx)){
  sample <- pheno$Sample_ID #pheno$Individual
  
  lmer.null  <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       pheno$DNAm_SV1 + pheno$DNAm_SV2 + pheno$DNAm_SV3 + 
                       pheno$DNAm_SmokingScore +
                       (1|sample), REML = F)
  
  lmer.model <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       pheno$DNAm_SV1 + pheno$DNAm_SV2 + pheno$DNAm_SV3 +  
                       pheno$DNAm_SmokingScore +
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
  
 # fdr.bh   <- round(p.adjust(res.anova$"Pr(>Chisq)"[2], method = "fdr"), 4) 
 # fdr.bon  <- round(p.adjust(res.anova$"Pr(>Chisq)"[2], method = "bonferroni"), 4)
    
  cbind(as.character(rownames(beta.mtrx)[cpg]), 
        signif(res.anova$"Pr(>Chisq)"[2], 4),
        round(res.anova$Chisq[2], 4),
        mean.veh, sd.veh, var.veh,
        mean.dex, sd.dex, var.dex,
        fc, var.all, mad)
}

stopImplicitCluster()

res.colnames <- cbind("PROBE_ID", "PVAL", "CHI_SQ",
                      "MEAN_VEH", "SD_VEH", "VAR_VEH",
                      "MEAN_DEX", "SD_DEX", "VAR_DEX",
                      "FC", "VAR_ALL", "MAD")
res          <- rbind(res.colnames, res)

write.table(res,
            file = lmer.res.out.fn, row.names = F, quote = F, sep = "\t", col.names = F, append = F)