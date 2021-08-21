# Prepare pheno

pheno.fn  <- "~/bio/datasets/pheno/pheno_full_for_kimono.csv"
pheno     <- read.csv2(pheno.fn, na.strings = "NA") 

dnam.svs <- read.csv2("~/bio/datasets/methylation/15_DNAm_sva/methyl_svs.csv")
summary(dnam.svs)

colnames(pheno) 
pheno <- left_join(pheno[, -c(77, 78, 79)], dnam.svs, by = "DNAm_ID" )
# pheno <- left_join(pheno, dnam.svs)
summary(pheno)

write.csv2(pheno, 
           pheno.fn, 
           quote = F, row.names = F)
