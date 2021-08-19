
sample.ids.fn <- "~/bio/code/mpip/dex-stim-human-array/data/sample_ids_for_analysis.csv"
sample.ids <- read.csv2(sample.ids.fn)

# Methylation data

src.data.dir  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/"
src.data.dir  <- "~/bio/datasets/methylation/10_final_qc_data/"
beta.mtrx.fn  <- paste0(src.data.dir, "dex_methyl_beta_combat_mtrx.rds")

beta.mtrx <- readRDS(beta.mtrx.fn)

beta.mtrx.new <- beta.mtrx[, colnames(beta.mtrx) %in% sample.ids$DNAm_ID]

saveRDS(beta.mtrx.new, 
        file = "~/bio/code/mpip/dex-stim-human-array/data/methylation/dex_methyl_beta_combat_mtrx.rds")


# Bio data

# Load full pheno tbl
pheno.full.fn <- paste0("~/bio/datasets/pheno/pheno_full_for_kimono.csv") 
pheno.full    <- read.csv2(pheno.full.fn) %>% setDT()
pheno         <- pheno.full[DNAm_ID %in% sample.ids$DNAm_ID]

write.csv2(pheno, 
           file = "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv", 
           row.names = F, quote = F)
