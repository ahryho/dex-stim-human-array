if (!require("MatrixEQTL")) 
  install.packages("MatrixEQTL")
library(MatrixEQTL)
library(data.table)

args <- commandArgs(T)

treatment    <- as.character(args[1]) #"veh"
eqtm.in.pre  <- as.character(args[2]) #"~/bio/datasets/eQTM/"
eqtm.res.pre <- as.character(args[3]) # paste0("~/bio/datasets/eQTM/result/")

treatment <- "dex"
type      <- "" # "_dnam_bcc" #""
mval <- "_beta" #"_mval"

# eqtm.in.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/test/"
# eqtm.res.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

eqtm.in.pre  <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
eqtm.res.pre <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

cpg.loc.fn  <- paste0(eqtm.in.pre, "cpg_locations.csv")
ensg.loc.fn <- paste0(eqtm.in.pre, "ensg_locations.csv")
snp.loc.fn  <- paste0(eqtm.in.pre, "snp_locations.csv")

# gex.layer.fn <- paste0(eqtm.in.pre, "gex_mtrx_", treatment, ".csv")

snp.layer.fn    <- paste0(eqtm.in.pre, "snp_mtrx.csv")
methyl.layer.fn <- paste0(eqtm.in.pre, "methyl", mval, "_mtrx_", treatment, ".csv")

if(treatment == "delta") 
  bio.layer.fn <- SlicedData$new() else
    bio.layer.fn <- paste0(eqtm.in.pre, "bio_mtrx_methyl_", treatment, type, ".csv")

eqtm.cis.result.fn <- paste0(eqtm.res.pre, "me-qtl_cis_result_", treatment, type, mval, ".csv")
eqtm.trans.result.fn <- paste0(eqtm.res.pre, "me-qtl_trans_result_", treatment, type, mval, ".csv")

# eqtm.cis.result.fn <- NULL
# eqtm.trans.result.fn <- NULL
# Load data

cpg.loc  <- fread(cpg.loc.fn)
# ensg.loc <- fread(ensg.loc.fn)
snp.loc  <- fread(snp.loc.fn)

# gex.layer <- fread(gex.layer.fn)
# snp.layer <- fread(snp.layer.fn)
# methyl.layer <- fread(methyl.layer.fn)
# bio.layer <- fread(bio.layer.fn)

# Check the colnames of layers are in the same order

# all(colnames(snp.layer)[-1] == colnames(methyl.layer)[-1])
# snp.layer[, colnames(snp.layer) %in% colnames(methyl.layer)]
# bio.layer[, colnames(bio.layer) %in% colnames(methyl.layer)]

RunMatrixEQTL <- function(snp.fn, gex.fn, bio.fn, cis.res.fn, trans.res.fn, cis.cutoff, trans.cutoff){
  
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
  snps                    <-  SlicedData$new()
  snps$fileDelimiter      <- ";"   # the TAB character
  snps$fileOmitCharacters <- "NA"  # denote missing values;
  snps$fileSkipRows       <- 1     # one row of column labels
  snps$fileSkipColumns    <- 1     # one column of row labels
  snps$fileSliceSize      <- 1e5   # read file in pieces of 2,000 rows
  snps$LoadFile(snp.fn)
  
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
  cvrt$fileSkipRows       <- 2    # one row of column labels + one row of treatment
  cvrt$fileSkipColumns    <- 1    # one column of row labels
  
  if(length(bio.fn) > 0){
    cvrt$LoadFile(bio.fn)  # read file if given
  }
  
  # 3. Run Matrix_eQTL
  me.all <- Matrix_eQTL_main(
    snps = snps1,
    gene = gene1,
    cvrt = cvrt1,
    useModel = useModel,
    errorCovariance = errorCovariance,
    snpspos = snp.loc, #cpg.loc,
    genepos = cpg.loc, #ensg.loc,
    output_file_name.cis = cis.res.fn,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    cisDist = cisDist,
    output_file_name = trans.res.fn,
    pvOutputThreshold = pvOutputThreshold_tra , #pvOutputThreshold_tra only cis
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    verbose = TRUE)
  
  return (me)
}

# Run matrixEQTL w
me.all <- RunMatrixEQTL(snp.fn = snp.layer.fn, # methyl.layer.fn, 
                        gex.fn = methyl.layer.fn,  # gex.layer.fn, 
                        bio.fn = bio.layer.fn, 
                        cis.res.fn = eqtm.cis.result.fn, 
                        trans.res.fn = eqtm.trans.result.fn, 
                        cis.cutoff = 5e-2, trans.cutoff = 0)

saveRDS(me.all, file =  paste0(eqtm.res.pre, "me-qtl_matrx_", treatment, type, mval, ".RDS"))
