

RunEnrichmentChromHMM <- function(own, background, chromhmm.states, fun, out.fn, nperm = 1000){
  
  no.cores <- detectCores() - 2
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  
  states.lst <- elementMetadata(chromhmm.states)[, "type"] %>% unique() %>% sort()
  
  enrich.perm.rslt <- foreach(i =  seq_along(states.lst), 
                                   .combine = rbind, 
                                   .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                         state <- states.lst[i]                                     
                                         public <- chromhmm.states[(elementMetadata(chromhmm.states)[, "type"]) == state, ] 
                                         fun(own = own, 
                                             background = background, 
                                             public = public, 
                                             nperm = nperm) 
                                       }
  
  stopImplicitCluster()
  
  enrich.perm.rslt <- cbind(enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                     state = states.lst)
  
  enrich.perm.rslt[["n_perm"]] <- nperm
  
  write.csv2(enrich.perm.rslt, 
             file = out.fn, 
             row.names = F, quote = F)
  
  return(enrich.perm.rslt)
}
