library(data.table)
library(dplyr)

src.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/"

all.eqtm.fn <- paste0(out.dir.pre, "eqtm_cis_all_fdr_005.csv")

all.df <- fread(all.eqtm.fn)

eqtms  <- list(delta = all.df[treatment == "delta", eQTM_ID], 
              dex = all.df[treatment == "dex", eQTM_ID], 
              veh = all.df[treatment == "veh", eQTM_ID])

only.dex.eqtms <- setdiff(setdiff(eqtms$dex, eqtms$veh), eqtms$delta)

only.dex.eqtms.df <- all.df[eQTM_ID %in% only.dex.eqtms]

fwrite(only.dex.eqtms.df,
       paste0(out.dir.pre, "eqtm_cis_unique_dex_fdr_005.csv"),
       quote = F, row.names = F, sep   = "\t")
