library(data.table)
library(dplyr)

source("~/bio/code/mpip/dex-stim-human-array/code/integrative/util.R")

src.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/"

eqtm.fn <- paste0(out.dir.pre, "eqtm_cis_all_fdr_005.csv")
eqtm.df <- fread(eqtm.fn)

eqtms  <- list(delta = eqtm.df[treatment == "delta", eQTM_ID], 
               dex = eqtm.df[treatment == "dex", eQTM_ID], 
               veh = eqtm.df[treatment == "veh", eQTM_ID])

only.dex.eqtms <- setdiff(setdiff(eqtms$dex, eqtms$veh), eqtms$delta)

only.dex.eqtms.df <- eqtm.df[eQTM_ID %in% only.dex.eqtms]

fwrite(only.dex.eqtms.df,
       paste0(out.dir.pre, "eqtm_cis_unique_dex_fdr_005.csv"),
       quote = F, row.names = F, sep   = "\t")
