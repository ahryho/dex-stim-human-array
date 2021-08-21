library(plyr)
library(dplyr)
library(limma)

# Differential methylation analysis

src.data.dir  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/"
beta.mtrx.fn  <- paste0(src.data.dir, "dex_methyl_beta_combat_mtrx.rds")
mval.mtrx.fn  <- paste0(src.data.dir, "dex_methyl_mval_mtrx_final.rds")
pheno.fn      <- paste0(src.data.dir, "dex_methyl_phenotype.rds")
cell.count.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/11_cell_types_estimation/dex_stim_array_human_cellcounts.csv"
dnam.age.fn   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/12_DNAm_age/dex_stim_array_human_meth_age_veh.csv"

dmr.data.dir <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMR/"

beta.mtrx <- readRDS(beta.mtrx.fn)
pheno     <- readRDS(pheno.fn)

# 1. Convert beta values to M-values

mval.mtrx <- t(apply(beta.mtrx, 1, function(x) log2(x / (1 - x))))
saveRDS(mval.mtrx, file = mval.mtrx.fn, compress = FALSE)

# 2. Prepare the pheno data

pheno     <- as.data.frame(pheno)
pheno.dma <- pheno %>% select(Individual = person, Sample_Name, Status = status, 
                              Sex = sex, Group = Sample_Group, Real_Age = age, BMI = bmi)

pheno.dma$Individual <- as.factor(pheno.dma$Individual)
pheno.dma$Sex        <- as.factor(pheno.dma$Sex)
pheno.dma$Group      <- as.factor(pheno.dma$Group)
pheno.dma$BMI        <- as.numeric(pheno.dma$BMI)

# 3. Add DNAm age as covariates

dnam.age.mtrx <- read.csv(dnam.age.fn, sep = ";", header = T )
dnam.age.mtrx <- as.data.frame(dnam.age.mtrx %>% select(Individual, 
                                                        DNAm_Age = PhenoAge))
pheno.dma     <- left_join(pheno.dma, dnam.age.mtrx, by = "Individual") 

pheno.dma[pheno$Sample_Name == "200705940062_R05C01", "DNAm_Age" ] <- 39

# 4. Add cell counts as covariates

cell.count.mtrx <- read.csv(cell.count.fn, sep = ";", header = T )
cell.count.mtrx <- cell.count.mtrx %>% mutate(Sample_Name = row.names(cell.count.mtrx))
pheno.dma       <- inner_join(pheno.dma, cell.count.mtrx, by = "Sample_Name") 

# 5. Making sure about samples in clinical and matrixes and their order

rownames(pheno.dma) <- pheno.dma$Sample_Name
 
table(colnames(mval.mtrx) %in% row.names(pheno.dma))
table(colnames(beta.mtrx) %in% row.names(pheno.dma))
#
all(row.names(pheno.dma) == colnames(mval.mtrx))
all(row.names(pheno.dma) == colnames(beta.mtrx))

# 6. Save results

write.table(pheno.dma, file = paste0(dmr.data.dir, "pheno.csv"), col.names = T, row.names = F, quote = F, sep = ";")

# 7. Extra
# sample.ids.fn <- "~/bio/code/mpip/dex-stim-human-array/data/sample_ids_for_analysis.csv"
# sample.ids <- read.csv2(sample.ids.fn)
# 
# src.data.dir  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/"
# src.data.dir  <- "~/bio/datasets/methylation/10_final_qc_data/"
# beta.mtrx.fn  <- paste0(src.data.dir, "dex_methyl_beta_combat_mtrx.rds")
# 
# beta.mtrx <- readRDS(beta.mtrx.fn)
# 
# beta.mtrx.new <- beta.mtrx[, colnames(beta.mtrx) %in% sample.ids$DNAm_ID]
# 
# saveRDS(beta.mtrx.new, 
#         file = "~/bio/code/mpip/dex-stim-human-array/data/methylation/dex_methyl_beta_combat_mtrx.rds")

