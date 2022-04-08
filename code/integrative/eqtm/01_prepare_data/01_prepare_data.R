# For initial data processing, please refer to the code in:
# "~/bio/code/mpip/dex-stim-human-array/code/integrative/meqtl/01_prepare_data/01_prepare_data.R"
library(data.table) 

GetDeltaMtrx <- function(mtrx.veh, mtrx.dex, id.col.idx = 1, out.fn){
  mtrx.delta <- mtrx.veh[,-1] - mtrx.dex[,-1] 
  mtrx.delta <- cbind(mtrx.veh[,1], mtrx.delta)
  
  fwrite(mtrx.delta, out.fn, quote = F, row.names = F, sep = ";")
  
  mtrx.delta
}

# Get GEX DeLTa mtrx from GEX DeX and GEX VeH residulas

lmer.res.out.fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/"

rslt.dir        <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"

gex.mtrx.veh    <- fread(paste0(lmer.res.out.fn, "gex_residuals/gex_svs_only/gex_residuals_veh.csv"))
gex.mtrx.dex    <- fread(paste0(lmer.res.out.fn, "gex_residuals/gex_svs_only/gex_residuals_dex.csv"))

all(rownames(gex.mtrx.veh) == rownames(gex.mtrx.dex))
order.idx  <- match(colnames(gex.mtrx.dex), colnames(gex.mtrx.veh))

all(colnames(gex.mtrx.veh) == colnames(gex.mtrx.dex))

gex.mtrx.delta <- GetDeltaMtrx(gex.mtrx.veh, gex.mtrx.dex, 1, out.fn = paste0(rslt.dir, "mtrx_delta_residuals/gex_mtrx_delta.csv"))

# Delta without cov and residuals
# 
# Delta GEX
gex.mtrx.veh   <- fread(paste0(rslt.dir, "gex_mtrx_veh.csv"))
gex.mtrx.dex   <- fread(paste0(rslt.dir, "gex_mtrx_dex.csv"), select = colnames(gex.mtrx.veh))

gex.mtrx.delta <- GetDeltaMtrx(gex.mtrx.veh, gex.mtrx.dex, 1, out.fn = paste0(rslt.dir, "gex_mtrx_delta.csv"))

# Delta DNAm
dnam.mtrx.veh   <- fread(paste0(rslt.dir, "methyl_beta_mtrx_veh.csv"))
dnam.mtrx.dex   <- fread(paste0(rslt.dir, "methyl_beta_mtrx_dex.csv"), select = colnames(dnam.mtrx.veh))

dnam.mtrx.delta <- GetDeltaMtrx(dnam.mtrx.veh, dnam.mtrx.dex, 1, out.fn = paste0(rslt.dir, "methyl_beta_mtrx_delta.csv"))

dnam_delta_check <- fread(paste0(rslt.dir, "methyl_beta_mtrx_delta_not_residuals.csv"))

# Delta DNAm residuals
dnam.mtrx.veh   <- fread(paste0(lmer.res.out.fn, "dnam_residuals/dnam_svs_only/dnam_residuals_veh.csv"))
dnam.mtrx.dex   <- fread(paste0(lmer.res.out.fn, "dnam_residuals/dnam_svs_only/dnam_residuals_dex.csv"))

all(rownames(dnam.mtrx.veh) == rownames(dnam.mtrx.dex))

dnam.mtrx.delta <- GetDeltaMtrx(dnam.mtrx.veh, dnam.mtrx.dex, 1, out.fn = paste0(rslt.dir, "mtrx_delta_residuals/methyl_beta_mtrx_delta.csv"))

dnam_delta_check <- fread(paste0(rslt.dir, "methyl_beta_mtrx_delta_not_residuals.csv"))
