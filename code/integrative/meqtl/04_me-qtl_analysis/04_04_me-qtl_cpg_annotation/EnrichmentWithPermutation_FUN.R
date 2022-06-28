
EnrichmentWithPermutationWithoutMAF <- function(own, background, public, nsample, nperm){
  
  # own <- meqtls.cpg.delta.coord.gr
  # background <- meqtls.cpg.veh.coord.gr
  # public <- anno.epic.gr[(elementMetadata(anno.epic.gr)[, "Relation_to_Island"]) == "N_Shelf", ]
  
  GenomeInfoDb::seqlevelsStyle(own) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(public) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(background) <- "UCSC"
  
  overlap     <- sum(overlapsAny(own, public))  # number of overlaps between the delta snps and public data (chromHMM_state or GWAS)
  # sum(distanceToNearest(own, public)@elementMetadata@listData$distance == 0)
  non_overlap <- length(own) - overlap
  
  overlap_bkgr      <- sum(overlapsAny(background, public))  # number of overlaps between baseline snps and public data (chromHMM_state or GWAS)
  non_overlap_bkgr  <- length(background) - overlap_bkgr
  
  conf_mtrx        <- matrix(c(overlap, overlap_bkgr, non_overlap, non_overlap_bkgr), 2, 2, byrow = TRUE)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  

  resampling <- lapply(1:nperm, function(x){
    sample_overlap     <- sum(IRanges::overlapsAny(sample(background, nsample), public))
    sample_non_overlap <- length(own) - sample_overlap
    conf_mtrx          <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
    fisher_test_rslt   <- fisher.test(conf_mtrx)
    c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate)
  }
  )

  p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] <= p_value) + 1) / (nperm + 1)
  p_value_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  
  # c(n_snps_overlap = overlap, or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp)
  
  return(c(n_overlap = overlap, or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp)) 
}

EnrichmentWithPermutationWithoutMAFnoPubData <- function(own, background, feature, nperm){
  
  GenomeInfoDb::seqlevelsStyle(own) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(background) <- "UCSC"
  
  nsample <- length(own)
  
  overlap      <- length(own[(elementMetadata(own)[, "annotation"]) == feature, ])  
  non_overlap  <- length(own) - overlap
  
  overlap_bkgr      <- length(background[(elementMetadata(background)[, "annotation"]) == feature, ])  
  non_overlap_bkgr  <- length(background) - overlap_bkgr
  
  conf_mtrx        <- matrix(c(overlap, overlap_bkgr, non_overlap, non_overlap_bkgr), 2, 2, byrow = TRUE)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  
  
  resampling <- lapply(1:nperm, function(x){
    sample_gr          <- sample(background, nsample)     
    sample_overlap     <- length(sample_gr[(elementMetadata(sample_gr)[, "annotation"]) == feature, ] )
    sample_non_overlap <- nsample - sample_overlap
    conf_mtrx          <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
    fisher_test_rslt   <- fisher.test(conf_mtrx)
    c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate)
  }
  )
  
  p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] <= p_value) + 1) / (nperm + 1)
  p_value_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  
  # c(n_snps_overlap = overlap, or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp)
  
  return(c(n_overlap = overlap, or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp)) 
}

EnrichmentWithPermutationGeneLocWithoutMAF <- function(own, background, public, nsample = NULL, nperm = 100){

  GenomeInfoDb::seqlevelsStyle(own) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(public) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(background) <- "UCSC"
  
  overlap      <- sum(overlapsAny(own, public))  
  overlap_bkgr <- sum(overlapsAny(background, public))
  
  if (is.null(nsample)) nsample <- length(own)
  
  conf_mtrx        <- matrix(c(overlap, length(own), overlap_bkgr, length(background)), 2, 2, byrow = F)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  
  
  resampling <- lapply(1:nperm, function(x){
    sample_overlap     <- sum(IRanges::overlapsAny(sample(background, nsample), public))
    conf_mtrx          <- matrix(c(overlap, length(own), sample_overlap, nsample), 2, 2, byrow = F)
    fisher_test_rslt   <- fisher.test(conf_mtrx)
    c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate)
  }
  )
  
  p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] <= p_value) + 1) / (nperm + 1)
  p_value_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  
  return(c(n_overlap = overlap, or = or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp)) 
}
