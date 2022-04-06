library(data.table)

require(lme4)
require(foreign)

library(parallel)
library(foreach)
library(doParallel)

# 1. Load data

treatment    <- as.character(args[1]) 
eqtm.in.pre  <- as.character(args[2]) 
eqtm.res.pre <- as.character(args[3]) 

eqtm.in.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/"
eqtm.res.pre <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/"

# gex.layer.fn    <- paste0(eqtm.in.pre, "gex_mtrx_", treatment, ".csv")
# methyl.layer.fn <- paste0(eqtm.in.pre, "methyl", mval, "_mtrx_", treatment, ".csv")
# 
# if(treatment != "delta") 
#     bio.layer.fn <- paste0(eqtm.in.pre, "bio_mtrx_gex_", treatment, type, ".csv")


gex.dex.layer.fn    <- paste0(eqtm.res.pre, "gex_residuals/gex_residuals_dex.csv")
gex.veh.layer.fn    <- paste0(eqtm.res.pre, "gex_residuals/gex_residuals_veh.csv")

dnam.dex.layer.fn    <- paste0(eqtm.res.pre, "dnam_residuals/dnam_residuals_dex.csv")
dnam.veh.layer.fn    <- paste0(eqtm.res.pre, "dnam_residuals/dnam_residuals_veh.csv")

# Prepare mtrx

pheno.fn <- paste0(eqtm.in.pre, "pheno/pheno_full_for_kimono.csv")
pheno    <- read.csv2(pheno.fn)#, na.strings = "#N/A") 
pheno    <- pheno[pheno$Include == 1, ]

# Prepare GEX mtrx
# 
samples.veh.ids <- pheno$RNA_ID[pheno$Group == "veh"]
samples.dex.ids <- pheno$RNA_ID[pheno$Group == "dex"]

gex.veh  <- fread(gex.veh.layer.fn)
colnames(gex.veh) <- c("ENSG_ID", samples.veh.ids)

gex.dex  <- fread(gex.dex.layer.fn)
colnames(gex.dex) <- c("ENSG_ID", samples.dex.ids)

gex.mtrx <- cbind(gex.veh, gex.dex[, -1])

# Prepare methyl mtrx
# 
samples.veh.ids <- pheno$DNAm_ID[pheno$Group == "veh"]
samples.dex.ids <- pheno$DNAm_ID[pheno$Group == "dex"]

dnam.veh  <- fread(dnam.veh.layer.fn)
colnames(dnam.veh) <- c("CpG_ID", samples.veh.ids)

dnam.dex  <- fread(dnam.dex.layer.fn)
colnames(dnam.dex) <- c("CpG_ID", samples.dex.ids)

dnam.mtrx <- cbind(dnam.veh, dnam.dex[, -1])

# Prepare pheno
# 
pheno$Sample_ID  <- as.factor(pheno$Sample_ID)

# Load dex only eqtms
# 
selected.df <- fread(paste0(eqtm.res.pre, "/eqtms/eqtm_cis_unique_dex_fdr_005.csv"))

# 2. Making sure about samples in pheno and and betas matrix in the same order

pheno   <- pheno[match(colnames(dnam.mtrx)[-1], pheno$DNAm_ID ),]
all(pheno$DNAm_ID == colnames(dnam.mtrx)[-1])
all(pheno$RNA_ID == colnames(gex.mtrx)[-1])

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# res <- foreach(cpg =  1:6, .combine = rbind, .packages = 'lme4') %dopar% {
res <- foreach(eqtm = 1:nrow(selected.df), .combine = rbind, .packages = 'lme4') %dopar% {
  
  sample  <- pheno$Sample_ID
  eqtm.id <- selected.df[eqtm]
  cpg    <- selected.df$CpG_ID[eqtm]
  gene   <- selected.df$ENSG_ID[eqtm]
  
  dnam <- t(as.matrix(dnam.mtrx[CpG_ID == cpg, -1]))
  gex  <- t(as.matrix(gex.mtrx[ENSG_ID == gene, -1]))
  
  lmer.null  <- lmer(gex ~ dnam + (1|sample), REML = F)
  lmer.model <- lmer(gex ~ dnam +  pheno$Group + (1|sample), REML = F)
  
  res.anova  <- anova(lmer.null, lmer.model)
  
  # Statistics
  # mean.veh <- round(mean(methyl.mtrx[cpg, samples.veh.ids]), 4)
  # mean.dex <- round(mean(methyl.mtrx[cpg, samples.dex.ids]), 4)
  # sd.veh   <- round(sd(methyl.mtrx[cpg, samples.veh.ids]), 4)
  # sd.dex   <- round(sd(methyl.mtrx[cpg, samples.dex.ids]), 4)
  # var.veh  <- round(var(methyl.mtrx[cpg, samples.veh.ids]), 4)
  # var.dex  <- round(var(methyl.mtrx[cpg, samples.dex.ids]), 4)
  
  # fc       <- round(mean.dex - mean.veh, 4)# round(log2(mean.dex) - log2(mean.veh), 4)
  # var.all  <- round(var(beta.mtrx[cpg, ]), 4)
  # mad      <- round(mad(beta.mtrx[cpg, ]), 4)
  
  cbind(selected.df$eQTM_ID[eqtm], 
        cpg, 
        gene, 
        p_value = signif(res.anova$"Pr(>Chisq)"[2], 4),
        chi_sq = round(res.anova$Chisq[2], 4))
}

stopImplicitCluster()

res[["fdr"]] <- p.adjust(res[, "p_value"], method = "fdr")
res <- cbind(res, fdr)

res.colnames <- c("eQTM_ID", "CpG_ID", "ENSG_ID", "PVAL", "CHI_SQ", "FDR")
res          <- rbind(res.colnames, res)

eqtm.cis.result.fn   <- paste0(eqtm.res.pre, "eqtms/eqtm_cis_result_unique_dex_fdr_005.csv")

write.table(res,
            file = eqtm.cis.result.fn, row.names = F, quote = F, sep = "\t", col.names = F, append = F)
