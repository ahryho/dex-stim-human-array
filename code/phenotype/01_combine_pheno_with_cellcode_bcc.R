library(dplyr)
library(data.table)

# Load full pheno tbl
# pheno.full.fn <- "~/bio/datasets/pheno/pheno_full_for_kimono.csv"
pheno.full.fn <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
pheno.full    <- read.csv2(pheno.full.fn)

# Load cell code cell bcc
bcc.cell.code.fn <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/cellcode.txt"
bcc.cell.code    <- read.table(bcc.cell.code.fn, sep = "\t", header = T)

# Inside cellcode
colnames(bcc.cell.code)
summary(bcc.cell.code)

bcc.cell.code <- bcc.cell.code[, c("RNA_ID", "Neutrophil", "Tcell", "Monocyte", "Bcell", "NKcell", "PlasmaCell", "DendriticCell")]
colnames(bcc.cell.code) <- c("RNA_ID", "Neutrophil_CellCODE", "Tcell_CellCODE", "Monocyte_CellCODE", 
                             "Bcell_CellCODE", "NKcell_CellCODE", "PlasmaCell_CellCODE", "DendriticCell_CellCODE")

pheno.full  <- left_join(pheno.full, bcc.cell.code, by = "RNA_ID")
colnames(pheno.full)

pheno.full <- dplyr::relocate(pheno.full, Include, .after = last_col())
colnames(pheno.full)

# Save 
write.csv2(pheno.full, 
           file = "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv", 
           row.names = F, quote = F)
