library(dplyr)
library(data.table)

# Load full pheno tbl
# pheno.full.fn <- "~/bio/datasets/pheno/pheno_full_for_kimono.csv"
pheno.full.fn <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
pheno.full    <- read.csv2(pheno.full.fn) %>% setDT()

# Load cell code cell bcc
bcc.fn <- "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/dex_stim_array_human_epidish_bcc.csv"
bcc    <- read.csv2(bcc.fn)

# Inside epidish bcc
colnames(bcc)
summary(bcc)

pheno.full  <- left_join(pheno.full, bcc, by = "DNAm_ID")
colnames(pheno.full)

pheno.full <- dplyr::relocate(pheno.full, Include, .after = last_col())
colnames(pheno.full)

pheno        <- pheno.full[Include == 1][Dex == 0]
bcc          <- pheno[, c("Cort_D1", "ACTH_D1", "Leukos_D1", "Gran_D1", "Mono_D1", "Lymph_D1")]
bcc          <- bcc %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
bcc.meth     <- pheno[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]
bcc.cellcode <- pheno[, c("Neutrophil_CellCODE", "Tcell_CellCODE", "Monocyte_CellCODE", "Bcell_CellCODE",
                          "NKcell_CellCODE", "PlasmaCell_CellCODE", "DendriticCell_CellCODE")]
bcc.epidish  <- pheno[, c("epidish_rpc.B", "epidish_rpc.NK", "epidish_rpc.CD4T", "epidish_rpc.CD8T", "epidish_rpc.Mono", "epidish_rpc.Neutro", "epidish_rpc.Eosino")]

sv.meth  <- pheno[, paste0("DNAm_SV", 1:3)]
snp.pcs  <- pheno[, paste0("PC", 1:2)]
sv.gex   <- pheno[, paste0("V", 1:5)] 

# Plot

tmp.df <- data.frame(#  bcc, bcc.cellcode,
  minfi = bcc.meth,
  bcc.epidish,
  sv.meth,
  snp.pcs,
  gex = sv.gex, 
  # Dex = as.numeric(pheno$Dex), 
  Status = as.numeric(pheno$Status), 
  Sex = as.numeric(pheno$Sex), BMI = pheno$BMI_D1, 
  SS  = pheno$DNAm_SmokingScore,
  Age = pheno$Age) %>% 
  na.omit() %>% 
  scale()

cor.mtrx <- cor(tmp.df)

library(corrplot)
corrplot(cor.mtrx, 
         method = "color", cl.pos = "r", cl.cex = 0.6,
         type = "upper",
         addgrid.col = "gray90",
         addCoef.col = "black", number.cex = 0.5, number.digits = 1, # Add coefficient of correlation
         tl.col = "red",  tl.cex = 0.6, tl.srt = 45, tl.offset = 0.5, tl.pos = "td", #Text label color and rotation
         sig.level = 0.01) 

# Save 
# write.csv2(pheno.full, 
#            file = "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv", 
#            row.names = F, quote = F)
