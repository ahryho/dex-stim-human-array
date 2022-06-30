library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/"

meqtl.veh.snp.gr   <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_veh_snps_with_maf_gr.rds"))
meqtl.delta.snp.gr <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_delta_snps_with_maf_gr.rds"))
background.all.gr  <- readRDS(paste0(out.dir.pre, "data/snps/imputed_qc/from_janine/qc/dex_geno_snps_with_maf_bins_gr.rds"))

# GWAS CD2
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds"))
# gwas.gw.gr <- gwas.gr[elementMetadata(gwas.gr)[, "p_value"] <= 5e-8, ] 

nperm <- 1000

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_CD_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.cd <- gwas.enrich.perm.rslt

# GWAS ADHD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_ADHD_Demontis_2019_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_ADHD_Demontis_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.adhd <- gwas.enrich.perm.rslt

# GWAS ASD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_ASD_Grove_2019_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_ASD_Grove_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.asd <- gwas.enrich.perm.rslt

# GWAS BPD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_BPD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_BP_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.bpd <- gwas.enrich.perm.rslt

# GWAS MDD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_MDD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_MDD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.mdd <- gwas.enrich.perm.rslt

# GWAS SCZ
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_SCZ_2018_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_SCZ_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.scz <- gwas.enrich.perm.rslt

# GWAS BMI
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/BMI_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_BMI_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.bmi <- gwas.enrich.perm.rslt

# GWAS IBD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/IBD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_IBD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.ibd <- gwas.enrich.perm.rslt

# GWAS PD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_PD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.pd <- gwas.enrich.perm.rslt

# filter GWAS CD datasets for genome-wide (gw) hits (p-value <= 5e-8) 
# 
# gwas.gw.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
#                                                       background = meqtl.veh.snp.gr, 
#                                                       public = gwas.gw.gr, 
#                                                       nperm = nperm)
# 
# gwas.gw.enrich.perm.rslt[["data"]] <- "GWAS_CD_2019_gw_p_5e-08"
# gwas.gw.enrich.perm.rslt[["n_perm"]] <- nperm
# gwas.gw.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gw.gr)

(rslt <- rbind(gwas.cd, gwas.adhd, gwas.asd, gwas.bpd, gwas.mdd, gwas.ibd, gwas.bmi, gwas.scz, gwas.pd))

# /05_me-qtl_enrichment/region_wise_independent_snps/
# 
write.csv2(rslt, 
           file = paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_GWAS_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# rslt <- read.csv2(file = paste0(out.dir.pre, "/05_me-qtl_enrichment/meqtl_snps_GWAS_CD_2019_enrichment_perm_", nperm, ".csv"))

# ChromHMM

chromhmm.blood.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_blood_states.Rds") 
chromhmm.brain.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_brain_states.Rds") 

states.lst <- elementMetadata(chromhmm.blood.states)[, "type"] %>% unique() %>% sort()

# chromHMM Blood

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 20

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
  state <- states.lst[i]                                     
  public <- chromhmm.blood.states[(elementMetadata(chromhmm.blood.states)[, "type"]) == state, ] 
  # rslt <- data.frame()
  # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
  EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                            background = meqtl.veh.snp.gr, 
                            public = public, 
                            nperm = nperm) 
}

 stopImplicitCluster()

chromhmm.enrich.perm.rslt <- cbind(chromhmm.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                   state = states.lst)

chromhmm.enrich.perm.rslt[["data"]] <- "Blood_and_T-cells"
chromhmm.enrich.perm.rslt[["n_perm"]] <- nperm

chromhmm.enrich.perm.rslt

# /05_me-qtl_enrichment/region_wise_independent_snps/

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "meqtl_snps_chromHMM_blood_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# chromHMM Brain

no.cores <- detectCores() - 3
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 20

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- states.lst[i]                                     
                                       public <- chromhmm.brain.states[(elementMetadata(chromhmm.brain.states)[, "name"]) == state, ] 
                                       # rslt <- data.frame()
                                       # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
                                       EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                                 background = meqtl.veh.snp.gr, 
                                                                 public = public, 
                                                                 nperm = nperm) 
                                     }

stopImplicitCluster()

chromhmm.enrich.perm.rslt <- cbind(chromhmm.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                   state = states.lst)

chromhmm.enrich.perm.rslt[["data"]] <- "Brain"
chromhmm.enrich.perm.rslt[["n_perm"]] <- nperm

chromhmm.enrich.perm.rslt

# /05_me-qtl_enrichment/region_wise_independent_snps/

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "meqtl_snps_chromHMM_brain_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# delta vs all

# chromHMM Blood

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 20

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- states.lst[i]                                     
                                       public <- chromhmm.blood.states[(elementMetadata(chromhmm.blood.states)[, "type"]) == state, ] 
                                       # rslt <- data.frame()
                                       # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
                                       EnrichmentWithPermutation(own = meqtl.veh.snp.gr, 
                                                                 background = background.all.gr, 
                                                                 public = public, 
                                                                 nperm = nperm) 
                                     }

stopImplicitCluster()

chromhmm.enrich.perm.rslt <- cbind(chromhmm.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                   state = states.lst)

chromhmm.enrich.perm.rslt[["data"]] <- "Blood_and_T-cells"
chromhmm.enrich.perm.rslt[["n_perm"]] <- nperm

chromhmm.enrich.perm.rslt

# /05_me-qtl_enrichment/region_wise_independent_snps/

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "meqtl_snps_chromHMM_blood_enrichment_veh_vs_all_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# chromHMM Brain

no.cores <- detectCores() - 3
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 20

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- states.lst[i]                                     
                                       public <- chromhmm.brain.states[(elementMetadata(chromhmm.brain.states)[, "name"]) == state, ] 
                                       # rslt <- data.frame()
                                       # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
                                       EnrichmentWithPermutation(own = meqtl.veh.snp.gr, 
                                                                 background = background.all.gr, 
                                                                 public = public, 
                                                                 nperm = nperm) 
                                     }

stopImplicitCluster()

chromhmm.enrich.perm.rslt <- cbind(chromhmm.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                   state = states.lst)

chromhmm.enrich.perm.rslt[["data"]] <- "Brain"
chromhmm.enrich.perm.rslt[["n_perm"]] <- nperm

chromhmm.enrich.perm.rslt

# /05_me-qtl_enrichment/region_wise_independent_snps/

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "meqtl_snps_chromHMM_brain_enrichment_veh_vs_all_perm_", nperm, ".csv"), 
           row.names = F, quote = F)
