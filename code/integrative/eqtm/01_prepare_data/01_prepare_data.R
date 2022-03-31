# For initial data processing, please refer to the code in:
# "~/bio/code/mpip/dex-stim-human-array/code/integrative/meqtl/01_prepare_data/01_prepare_data.R"
 

# Get GEX DeLTa mtrx from GEX DeX and GEX VeH residulas

lmer.res.out.fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/gex_residuals/gex_residuals/"

gex.mtrx.veh    <- fread(paste0(lmer.res.out.fn, "gex_residuals_veh.csv"))
gex.mtrx.dex    <- fread(paste0(lmer.res.out.fn, "gex_residuals_dex.csv"), select = colnames(gex.mtrx.veh))

# methyl.mtrx.veh <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals_veh.csv")
# methyl.mtrx.dex <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals_dex.csv")

all(rownames(gex.mtrx.veh) == rownames(gex.mtrx.dex))
order.idx  <- match(colnames(gex.mtrx.dex), colnames(gex.mtrx.veh))

# methyl.mtrx.dex <- methyl.mtrx.dex[, colnames(methyl.mtrx.veh)]
all(colnames(gex.mtrx.veh) == colnames(gex.mtrx.dex))

gex.mtrx.delta <- gex.mtrx.veh[,-1] - gex.mtrx.dex[,-1] 
gex.mtrx.delta <- cbind(gex.mtrx.veh[,1], gex.mtrx.delta)

fwrite(gex.mtrx.delta, 
       paste0(lmer.res.out.fn, "gex_mtrx_delta.csv"),
       quote = F, row.names = F, sep = ";")
