# Cell-type enrichment
# 
# Load packages
#

library(TCA)
library("ggplot2")
library("ggpubr")
library("pracma")
library("matrixStats")
library(data.table)
library(dplyr) 

library(parallel)

# Load data
# 
setwd("~/bio/code/mpip/dex-stim-human-array/")
setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/")

# Beta mtrx
# 
methyl.mtrx.veh <- fread("data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv") 

# Pheno 
# 
pheno.full            <- fread("data/pheno/pheno_full_for_kimono.csv", dec = c(",")) 
bcc.epidish.rpc.salas <- read.csv2( file = "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv")
colnames(bcc.epidish.rpc.salas)[2:13] <- paste0("salas.", colnames(bcc.epidish.rpc.salas)[2:13])

pheno                 <- pheno.full[Include == 1][Dex == 0]
pheno                 <- left_join(pheno, bcc.epidish.rpc.salas)

cols   <- colnames(bcc.epidish.rpc.salas)[2:13]
bcc.df <- data.frame(pheno[, ..cols], row.names = pheno$DNA_ID) # W

cov.lst <- c(
  "Sex", "Status", "Age", "BMI_D1", "DNAm_SmokingScore",
  "PC1", "PC2")

cov.df <- data.frame(pheno[, ..cov.lst], row.names = pheno$DNA_ID)

methyl.mtrx <- data.frame(methyl.mtrx.veh, row.names = methyl.mtrx.veh$CpG_ID) %>% select(-CpG_ID)

#####
# # Group CpGs
# #
# out.dir.pre <- "output/data/integrative/matrixEQTL/meqtls/primary_meqtl_results/"
# 
# meqtl.parallel.fc.grp.veh.dex.gr <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_veh_dex_gr.rds"))
# meqtl.parallel.fc.grp.delta.gr   <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_delta_gr.rds"))
# meqtl.opposite.fc.gr             <- readRDS(file =  paste0(out.dir.pre, "meqtl_overlaps_meqtl_level_cpg_anno_opposite_fc_grp_gr.rds"))
# 
# # Prepare data
# # 
# parallel.fc.grp.veh.dex.cpg <- list(CpG_ID = names(meqtl.parallel.fc.grp.veh.dex.gr@anno))
# parallel.fc.grp.delta.cpg   <- list(CpG_ID = names(meqtl.parallel.fc.grp.delta.gr@anno))
# opposite.fc.grp.cpg         <- list(CpG_ID = names(meqtl.opposite.fc.gr@anno))

# Oposite 
# 
methyl.mtrx <- methyl.mtrx.veh[CpG_ID %in% opposite.fc.grp.cpg$CpG_ID, ] # X
methyl.mtrx <- data.frame(methyl.mtrx, row.names = methyl.mtrx$CpG_ID) %>% select(-CpG_ID)

# Run TCA
# ####

tca.mdl.veh <- tca(
  X = methyl.mtrx,
  W = bcc.df,
  C1 = cov.df,
  refit_W = T,
  parallel = T,
  num_cores = detectCores() - 1)

treatment <- "veh"
out.fn <- paste0("output/data/integrative/cell_type_enrichment/dnam_opposite_fc_cell_type_enrichment_", treatment, ".RDS")
saveRDS(tca.mdl.veh, file = out.fn)

mus <- tca.mdl.veh$mus_hat
pval <- tca.mdl.veh$gammas_hat_pvals

tca.mdl.veh$W

tca.mdl.hannum.pvals.joint <- tca.mdl.hannum$gammas_hat_pvals.joint[,"Age"]
# Extract p-values for each cell type for its association with age, under a marginal conditional test
tca.mdl.hannum.pvals.marg_cond <- 
  tca.mdl.hannum$gammas_hat_pvals[,paste(colnames(bcc.df), ".Age",sep="")]  
# qq-plots - for the p-values of the joint test, and for the marginal conditional p-values for CD4 cells
plot_qq <- function(pvals, labels, ggarrange.nrow = 1, ggarrange.ncol = 1, alpha = 0.05, experiment_wide_line = TRUE){
  significance_th <- list(alpha/length(pvals[[1]]))
  if(length(pvals)-1) significance_th[[2]] <- alpha/(length(pvals)*length(pvals[[1]]))
  qqplots <- lapply(1:length(pvals), function(p){
    df <- data.frame(pvals.obs = -log10(sort(pvals[[p]])), pvals.exp = -log10(sort((1:length(pvals[[1]]))/length(pvals[[1]]))));
    qqplot <- ggplot(df, aes(x = pvals.exp, y = pvals.obs)) +
      stat_binhex(geom = "point") + #, bins=10, size=1) +
      geom_abline() +
      ggtitle(labels[p]) +
      xlab(expression(Expected~-log[10](P))) + ylab(expression(Observed~-log[10](P))) +
      theme_bw() +
      guides(fill="none") +
      geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)
    if (length(significance_th)-1 & experiment_wide_line) qqplot <- qqplot + geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=0.5)
    return(qqplot)
  })
  ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
}

plot_qq(list(tca.mdl.hannum.pvals.joint, tca.mdl.hannum.pvals.marg_cond[,"salas.CD8nv.Age"]),
        labels = c("Joint test with age", "CD4 marginal conditional test with age"),
        ggarrange.nrow = 1,
        ggarrange.ncol = 2,
        experiment_wide_line = FALSE)

