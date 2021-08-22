library(dplyr)
library(data.table)

# Load full pheno tbl
# pheno.full.fn <- paste0("~/bio/datasets/pheno/pheno_full_for_kimono.csv") 
pheno.full.fn <- paste0("~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv") 
pheno.full    <- read.csv2(pheno.full.fn)

colnames(pheno.full)

# pheno.full <- pheno.full[ , -which(names(pheno.full) %in% c("PC1", "PC2"))]
# colnames(pheno.full)

# Include PC1 and PC2 from genotype
eigenvec.fn  <- "~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_darina/ld_200_100/pca/mds/dex_geno_imputed_ld_mds.mds"
eigenvec     <- fread(eigenvec.fn)
pheno.full   <- left_join(pheno.full, eigenvec[,.(FID, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10)], by = c("DNA_ID" = "FID"))
colnames(pheno.full)[145:154] <- paste0("P", colnames(pheno.full)[145:154])

# Include into the further analysis only samples which are present in all data type
pheno.full[["Include"]] <- ifelse(!is.na(pheno.full$DNAm_ID), "1", "0")

# Mark genotype outliers

outlier.list <- c("MPIPSYKL_004888", "MPIPSYKL_007801", "MPIPSYKL_008512", "MPIPSYKL_008648")
pheno.full[pheno.full$DNA_ID %in% outlier.list, "Include"] <- 0
table(pheno.full$Include)


# Save 
write.csv2(pheno.full, 
           file = "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv", 
           row.names = F, quote = F)
