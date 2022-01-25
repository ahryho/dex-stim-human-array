library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

meqtl.veh.snp.gr   <- readRDS(paste0(out.dir.pre, "meqtl_veh_snps_with_maf_gr.rds"))
meqtl.delta.snp.gr <- readRDS(paste0(out.dir.pre, "meqtl_delta_snps_with_maf_gr.rds"))

# GWAS
gwas.gr    <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds")
gwas.gw.gr <- gwas.gr[elementMetadata(gwas.gr)[, "p_value"] <= 5e-8, ] 

nperm <- 100

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_CD_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm

# filter GWAS datasets for genome-wide (gw) hits (p-value <= 5e-8) 

gwas.gw.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                      background = meqtl.veh.snp.gr, 
                                                      public = gwas.gw.gr, 
                                                      nperm = nperm)

gwas.gw.enrich.perm.rslt[["data"]] <- "GWAS_CD_2019_gw_p_5e-08"
gwas.gw.enrich.perm.rslt[["n_perm"]] <- nperm

(rslt <- rbind(gwas.enrich.perm.rslt, gwas.gw.enrich.perm.rslt))

write.csv2(rslt, 
           file = paste0(out.dir.pre, "/05_me-qtl_enrichment/meqtl_snps_GWAS_CD_2019_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# rslt <- read.csv2(file = paste0(out.dir.pre, "/05_me-qtl_enrichment/meqtl_snps_GWAS_CD_2019_enrichment_perm_", nperm, ".csv"))

# ChromHMM

chromhmm.blood.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_blood_states.Rds") 
chromhmm.brain.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_brain_states.Rds") 

states.lst <- elementMetadata(chromhmm.blood.states)[, "type"] %>% unique() %>% sort()

# chromHMM Blood

no.cores <- detectCores() - 3
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 1

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

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "/05_me-qtl_enrichment/meqtl_snps_chromHMM_blood_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# chromHMM Brain

no.cores <- detectCores() - 3
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm <- 1

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

write.csv2(chromhmm.enrich.perm.rslt, 
           file = paste0(out.dir.pre, "/05_me-qtl_enrichment/meqtl_snps_chromHMM_brain_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)