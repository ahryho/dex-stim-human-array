if (!require("MatrixEQTL")) 
  install.packages("MatrixEQTL")
library(MatrixEQTL)
library(data.table)

args <- commandArgs(T)

treatment    <- as.character(args[1]) #"veh"
eqtm.in.pre  <- as.character(args[2]) #"~/bio/datasets/eQTM/"
eqtm.res.pre <- as.character(args[3]) # paste0("~/bio/datasets/eQTM/result/")

# treatment <- "dex"
type      <- "" # "_dnam_bcc" #""
mval <- "_beta" #"_mval"

# eqtm.in.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
# eqtm.res.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtm/"
# 
# eqtm.in.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"
# eqtm.res.pre <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/"

cpg.loc.fn  <- paste0(eqtm.in.pre, "cpg_locations.csv")
ensg.loc.fn <- paste0(eqtm.in.pre, "ensg_locations.csv")

if(treatment == "delta"){
	gex.layer.fn    <- paste0(eqtm.in.pre, "mtrx_delta_residuals/gex_mtrx_", treatment, ".csv")
	methyl.layer.fn <- paste0(eqtm.in.pre, "mtrx_delta_residuals/methyl", mval, "_mtrx_", treatment, ".csv")}

# if(treatment == "delta") 
#  bio.layer.fn <- SlicedData$new() else
    bio.layer.fn <- paste0(eqtm.in.pre, "bio_mtrx_methyl_gex_", treatment, type, ".csv")

eqtm.cis.result.fn   <- paste0(eqtm.res.pre, "eqtm_cis_result_", treatment, type, mval, ".csv")
eqtm.trans.result.fn <- paste0(eqtm.res.pre, "eqtm_trans_result_", treatment, type, mval, ".csv")

# Load data

cpg.loc  <- fread(cpg.loc.fn)[, .(CpG_ID, chr, pos)]
ensg.loc <- fread(ensg.loc.fn)

# gex.layer    <- fread(gex.layer.fn) 
# methyl.layer <- fread(methyl.layer.fn)
# bio.layer    <- fread(bio.layer.fn)

# Check the colnames of layers are in the same order

# all(colnames(gex.layer)[-1] == colnames(methyl.layer)[-1])
# all(colnames(bio.layer)[-1] == colnames(methyl.layer)[-1])
# all(colnames(bio.layer)[-1] == colnames(gex.layer)[-1])

RunMatrixEQTL <- function(methyl.fn, gex.fn, bio.fn, cis.res.fn, trans.res.fn, cis.cutoff, trans.cutoff){
  
  # 1. Set up general parameters
  useModel              <- modelLINEAR #other models: modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  # We need to set up the p-value cutoff and the error model (in this case assuming independent errors)
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis <- cis.cutoff  # cis-eQTLS cutoff: 0.05 = 5e-2
  pvOutputThreshold_tra <- trans.cutoff  # trans-eQTLs cutoff: 0.01 = 1e-2 (0 means no trans-eQTLs)
  errorCovariance       <- numeric()
  cisDist               <- 1e6  # # Distance for local (cis) gene-SNP pairs: cis window of 1Mb = 1e6, 100kb = 1e5
  
  # 2. Data set up
  
  # SNP Data
  cpgs                    <-  SlicedData$new()
  cpgs$fileDelimiter      <- ";"   # the TAB character
  cpgs$fileOmitCharacters <- "NA"  # denote missing values;
  cpgs$fileSkipRows       <- 1     # one row of column labels
  cpgs$fileSkipColumns    <- 1     # one column of row labels
  cpgs$fileSliceSize      <- 1e5   # read file in pieces of 2,000 rows
  cpgs$LoadFile(methyl.fn)
  
  # GEX Data
  gene                    <- SlicedData$new()
  gene$fileDelimiter      <-  ";"  # the TAB character
  gene$fileOmitCharacters <- "NA"  # denote missing values;
  gene$fileSkipRows       <- 1     # one row of column labels
  gene$fileSkipColumns    <- 1     # one column of row labels
  gene$fileSliceSize      <- 1e5   # read file in pieces of 100,000 rows
  gene$LoadFile(gex.fn)
  
  # Biological data
  cvrt                    <- SlicedData$new()
  cvrt$fileDelimiter      <- ";"  # the TAB character
  cvrt$fileOmitCharacters <- "NA" # denote missing values;
  cvrt$fileSkipRows       <- 1    # one row of column labels
  cvrt$fileSkipColumns    <- 1    # one column of row labels
  
  if(length(bio.fn) > 0){
    cvrt$LoadFile(bio.fn)  # read file if given
  }
  
  # 3. Run Matrix_eQTL
  me.all <- Matrix_eQTL_main(
    snps = cpgs,
    gene = gene,
    cvrt = cvrt,
    useModel = useModel,
    errorCovariance = errorCovariance,
    snpspos = cpg.loc,
    genepos = ensg.loc,
    output_file_name.cis = cis.res.fn,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    cisDist = cisDist,
    output_file_name = trans.res.fn,
    pvOutputThreshold = pvOutputThreshold_tra , #pvOutputThreshold_tra only cis
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    verbose = TRUE)
  
  return (me.all)
}

# Run matrixEQTL w
me.all <- RunMatrixEQTL(methyl.fn = methyl.layer.fn, 
                        gex.fn = gex.layer.fn, 
                        bio.fn = bio.layer.fn, 
                        cis.res.fn = eqtm.cis.result.fn, 
                        trans.res.fn = eqtm.trans.result.fn, 
                        cis.cutoff = 5e-2, trans.cutoff = 0)

saveRDS(me.all, file =  paste0(eqtm.res.pre, "eqtm_matrx_", treatment, type, mval, ".RDS"))
