# Function to extract CpGs

GetMethylSubset <- function(treatment, public.data, is.out = T){
  beta.mtrx     <- LoadMethylBeta(treatment)
  beta.sub.mtrx <- beta.mtrx[CpG_ID %in% public.data$CpG_ID] #442
  
  if(is.out)
    fwrite(beta.sub.mtrx, 
           paste0("output/data/public_data/Provencal2019_Hipocampal_CpGs/mpip_methyl_beta_mtrx_", treatment, ".csv"),
           quote = F, row.names = F, sep = ";")
  
  return(beta.sub.mtrx)
}

# Function to extract cell type-specific analysis results

GetFDRSubset <- function(treatment, public.data){
  pval.df <- fread(paste0("output/data/integrative/cell_type_enrichment/dnam_cell_type_enrichment_pvals_", treatment, ".csv"))
  
  pval.sub.df <- pval.df[CpG_ID %in% public.data$CpG_ID, 1:13]
  fdr.bcc.df <- matrix(p.adjust(as.vector(as.matrix(pval.sub.df[, 2:13])), method ='fdr'), ncol = 12) %>% data.frame()
  fdr.bcc.df <- cbind(pval.sub.df$CpG_ID, fdr.bcc.df)
  colnames(fdr.bcc.df) <- colnames(pval.sub.df)
  
  fwrite(fdr.bcc.df, 
         paste0("output/data/public_data/Provencal2019_Hipocampal_CpGs/mpip_dnam_cell_type_enrichment_fdr_", treatment, ".csv"),
         quote = F, row.names = F, sep = ";")
  
  return(fdr.bcc.df)
}

GetEpiStressScore <- function(public.data, beta.mtrx, treatment){
  public.data <- public.data[public.data$CpG_ID %in%  beta.mtrx$CpG_ID, ]
  public.data <- public.data[match(beta.mtrx$CpG_ID, public.data$CpG_ID,), ]
  
  epistress.score <- apply(beta.mtrx[, -1], 2, function(sample){
    sum(sample * public.data$Weights)
  })
  
  epistress.score.df <- data.frame( DNA_ID = colnames(beta.dex.sub.mtrx[, -1]), EpiStressScore = epistress.score)
  
  fwrite(epistress.score.df, 
         paste0("output/data/public_data/Provencal2019_Hipocampal_CpGs/mpip_epistress_score_", treatment, ".csv"),
         quote = F, row.names = F, sep = ";")
  
}
