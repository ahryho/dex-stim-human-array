library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/"

# Relation to island

meqtls.cpg.delta.coord.gr <- readRDS(paste0(out.dir.pre, "meqtl_delta_cpgs_gr.rds"))
meqtls.cpg.veh.coord.gr   <- readRDS(paste0(out.dir.pre, "meqtl_veh_cpgs_gr.rds"))
anno.epic.gr              <- readRDS(paste0(out.dir.pre, "anno_epic_gr.rds"))    

features.lst <- elementMetadata(meqtls.cpg.delta.coord.gr)[, "Relation_to_Island"] %>% unique() %>% sort()

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 10
nsample <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- features.lst[i]                                     
                                       public <- anno.epic.gr[(elementMetadata(anno.epic.gr)[, "Relation_to_Island"]) == state, ] 
                                       EnrichmentWithPermutationWithoutMAF3(own = meqtls.cpg.delta.coord.gr,
                                                                           background = meqtls.cpg.veh.coord.gr,
                                                                           public = public,
                                                                           nsample = nsample,
                                                                           nperm = nperm)
                                     }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Faetute = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

# delta.meqtl.snp.anno.rds <- readRDS(paste0(out.dir.pre, "meqtls_cpg_annotated_withChIPseeker_delta.rds"))
# veh.meqtl.snp.anno.rds   <- readRDS(paste0(out.dir.pre, "meqtls_cpg_annotated_withChIPseeker_veh.rds"))
# 
# features.lst <- elementMetadata(delta.meqtl.snp.anno.rds@anno)[, "annotation"] %>% unique() %>% sort()

### chromHMM Blood
###

# meqtl.veh.gr   <- readRDS(paste0(out.dir.pre, "meqtl_veh_cpgs_gr.rds"))
# meqtl.delta.gr <- readRDS(paste0(out.dir.pre, "meqtl_delta_cpgs_gr.rds"))

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

# genome locations

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
