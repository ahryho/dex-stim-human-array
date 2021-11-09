
library(dplyr)
library(data.table)

# Set up parameters

snp.src.dir <- "~/bio/code/mpip/dex-stim-human-array/data/snps/final_imputed_qc_snps/mtrx/"
# snp.src.dir <- "/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/"

src.pheno.data.pre <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/"
# src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"

pheno.fn        <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")

# kimono.data.pre <- "~/bio/datasets/kimono/"
# kimono.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"
# Load data

x <- fread(paste0(snp.src.dir, "plinkgen_X.txt"))
y <- fread(paste0(snp.src.dir, "dex_geno_imputed.sample"), header = F)
z <- fread(paste0(snp.src.dir, "snps_positions.txt"), header = F)

pheno <- fread(pheno.fn, na.strings = c('#N/A', ''))
pheno <- pheno[Include == 1]

# Substract one sub-matrix from another: first number of each patient - second number 

#  1 indicates homozygote for reference
# -1 indicates heterozygote
#  0 inidicates homozygote for alternative

numsam <- (nrow(y) - 2) * 3
new_x  <- x[, seq(1, numsam, by = 3), with = F] - x[, seq(1, numsam, by = 3) + 1, with = F]
new_x[1:5]

# Code the numbers so that 0: homo-reference; 1: hetero; 2: homo-alternative
new_x[new_x == 1]  <-"ref"
new_x[new_x == -1] <-"hetero"
new_x[new_x == 0]  <-"alt"

new_x[new_x == "ref"]    <- 0
new_x[new_x == "hetero"] <- 1
new_x[new_x == "alt"]    <- 2

z$ID <- rownames(z)
colnames(z)[1] <- "SNP"
z$SNP           <- as.character(z$SNP)

new_x_1 <- cbind(z$ID, new_x)

# Add DNA_ID to the mat
colnames(new_x_1) <- c("sample", y$V1[3:length(y$V1)])

# Keep only columns (samples) that are in the rest of the data

# dna.ids   <- unique(pheno[!is.na(pheno$DNAm_ID), "DNA_ID"])
dna.ids   <- unique(pheno[Include == 1, "DNA_ID"])
whichcols <- c("sample", colnames(new_x_1)[colnames(new_x_1) %in% dna.ids$DNA_ID])

# extract them
new_x_2 <- new_x_1[,..whichcols, with = F] 

# as integer
changeCols <- whichcols[-1]
new_x_2[,(changeCols):= lapply(.SD, as.integer), .SDcols = changeCols]
str(new_x_2)

new_x_2 <- cbind(z$SNP, new_x_2)
new_x_2 <- new_x_2[, -c("sample")]

fwrite(new_x_2, 
       "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/snp_mtrx.csv",
       quote = F, row.names = F, sep = ";")
