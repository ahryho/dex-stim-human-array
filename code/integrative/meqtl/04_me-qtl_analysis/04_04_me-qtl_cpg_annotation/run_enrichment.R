library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/"

meqtl.veh.gr   <- readRDS(paste0(out.dir.pre, "meqtl_veh_cpgs_gr.rds"))
meqtl.delta.gr <- readRDS(paste0(out.dir.pre, "meqtl_delta_cpgs_gr.rds"))

# ChromHMM

chromhmm.blood.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_blood_states.Rds") 
chromhmm.brain.states <- readRDS("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/chromHMM_brain_states.Rds") 

states.lst <- elementMetadata(chromhmm.blood.states)[, "type"] %>% unique() %>% sort()

# chromHMM Blood

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 100
nsample <- 1000

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
  state <- states.lst[i]                                     
  public <- chromhmm.blood.states[(elementMetadata(chromhmm.blood.states)[, "type"]) == state, ] 
  # rslt <- data.frame()
  # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
  EnrichmentWithPermutationWithoutMAF(own = meqtl.delta.gr, 
                                      background = meqtl.veh.gr, 
                                      public = public, 
                                      nsample = nsample,
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
           file = paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/", clumping.approach, "_independent_snps/", "meqtl_cpgs_chromHMM_blood_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)

# chromHMM Brain

no.cores <- detectCores() - 3
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 100
nsample <- 1000

chromhmm.enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                     .combine = rbind, 
                                     .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- states.lst[i]                                     
                                       public <- chromhmm.brain.states[(elementMetadata(chromhmm.brain.states)[, "name"]) == state, ] 
                                       # rslt <- data.frame()
                                       # rslt[c("n_snps_overlap", "or", "or_perm", "p_val", "p_val_perm", "p_val_emp"), ] <- 
                                       EnrichmentWithPermutationWithoutMAF(own = meqtl.delta.gr, 
                                                                           background = meqtl.veh.gr, 
                                                                           public = public, 
                                                                           nsample = nsample, 
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
           file = paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/", clumping.approach, "_independent_snps/", "meqtl_cpgs_chromHMM_brain_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)
