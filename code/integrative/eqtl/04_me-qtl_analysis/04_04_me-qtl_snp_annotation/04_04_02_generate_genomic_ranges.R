library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(IRanges)
library(arules)
library(reshape2)
library(stringr)

# function to generate an eQTL GRanges object (with MAF bins)
# @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
# @param mdata Merging data (e.g. Background SNPs)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObjecteQTL <- function(input, mdata, ofile){
  # merge with MAF info
  input <- merge(input, mdata, by = "SNP")
  # generate GRanges object
  input_gr <- GenomicRanges::GRanges(seqnames = input$CHR, 
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)), 
                                                               end = as.numeric(as.character(input$POS))),
                                     snp_id = input$SNP, 
                                     bin = input$bin)
  # remove duplicates
  input_gr <- input_gr[!duplicated(input_gr)] 
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "UCSC"
  # save GRanges object to an rds file
  saveRDS(input_gr, file =  ofile)
  # return finalized GRanges object
  return(input_gr)
}

# function to generate a GRanges object (for all datasets with Chromosome, base pair position and rs SNP ID information)
# @param input Input data of which an GRanges object should be made of (needs CHR, BP, P and SNP as column names)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObject <- function(input, ofile){
  # filter for p-value (P) < 0.05
  input <- input[input$P < 0.05, ]
  # remove all rows including NA in any column
  input <- input[complete.cases(input), ]
  # generate GRanges object
  input_gr <- GenomicRanges::GRanges(seqnames = input$CHR,
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)),
                                                               end = as.numeric(as.character(input$POS))),
                                     snp_id = input$SNP,
                                     p_value = input$P)  # generate GRanges object only with chr, position, rs SNP ID and p-value
  # remove duplicates
  input_gr <- input_gr[!duplicated(input_gr)] 
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "UCSC"
  print(length(input_gr)) # print length
  # save GRanges object to an rds file
  saveRDS(input_gr, file = ofile)
  return(input_gr)
}


out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/region_wise_independent_snps/`"

# background snps (baseline + delta )
# Read MAF values for each SNP
snps.dir <- "~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_janine/qc/"
mafs     <- fread(paste0(snps.dir, "dex_geno_maf_for_filtered_samples.afreq"), data.table = F)
colnames(mafs) <- c("CHR", "SNP", "REF", "ALT", "MAF", "OBS_CT")

# generate bins of MAF in 0.05 steps
boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
mafs$bin <- arules::discretize(mafs$ALT_FREQS, method = "fixed", breaks = boundaries, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# load background SNPs
snps <- fread(paste0(snps.dir, "Dex_genoData_SNPs.bim"), data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos

background.all.gr <- GenerateGrangesObjecteQTL(input = snps, 
                                           mdata = mafs, 
                                           ofile = paste0(snps.dir, "dex_geno_snps_with_maf_bins_gr.rds"))

#############

background.all.gr <- readRDS( paste0(snps.dir, "dex_geno_snps_with_maf_bins_gr.rds"))

ind.meqtl.veh.df   <- fread(meqtl.veh.fn, col.names = col.names)
ind.meqtl.dex.df   <- fread(meqtl.dex.fn, col.names = col.names)
ind.meqtl.delta.df <- fread(meqtl.delta.fn, col.names = col.names)

# snps.bkgr <- snps[snps$SNP %in% ind.meqtl.delta.df$SNP,]
# 
# # Generate background SNP data as GRanges object
# background.gr <- GenerateGrangesObjecteQTL(input = snps.bkgr, 
#                                            mdata = mafs, 
#                                            ofile = paste0(out.dir.pre, "meqtl_delta_snps_with_maf_gr.rds"))
# table(seqnames(background.gr))  # no chr X/Y/M (remove in dataset if included there)
# length(background.gr$snp_id) 

# Generate merging file for eQTL data from background SNPs 
background <- as.data.frame(background.all.gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")

meqtl.delta.snp.gr <- GenerateGrangesObjecteQTL(input = ind.meqtl.delta.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "meqtl_delta_snps_with_maf_gr.rds"))

meqtl.veh.snp.gr   <- GenerateGrangesObjecteQTL(input = ind.meqtl.veh.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "meqtl_veh_snps_with_maf_gr.rds"))

meqtl.dex.snp.gr   <- GenerateGrangesObjecteQTL(input = ind.meqtl.dex.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "meqtl_dex_snps_with_maf_gr.rds"))
# GWAS Summary Statistics

gwas.cluster.out.dir <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/"

# gwas.cross.dis <- read.table(text = gsub(" ", "\t", readLines("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_10k_oct2019_v2.txt.daner.txt")), header = T)
# gwas.cross.dis.10k <- fread("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_10k_oct2019_v2.txt.daner.txt", 
#                             header = T, stringsAsFactors = F, 
#                             select = c("CHROM", "POS", "ID", "PVAL"))
# colnames(gwas.cross.dis.10k) <- c("CHR", "POS", "SNP", "P")

# Cross Disorders

gwas.cross.dis <- fread("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt", 
                        header = T, stringsAsFactors = F, 
                        select = c("CHROM", "POS", "ID", "PVAL"))
colnames(gwas.cross.dis) <- c("CHR", "POS", "SNP", "P")

gwas.gr <- GenerateGrangesObject(gwas.cross.dis, 
                                 ofile = "~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds")

# SCZ 

fn <- "/binder/common/public_data/PGC/SCZ_2018/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_SCZ_2018_GR_p005.rds"))

# ADHD
 
fn <- "/binder/common/public_data/PGC/ADHD/Demontis_2019/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_ADHD_Demontis_2019_GR_p005.rds"))

# ASD

fn <- "/binder/common/public_data/PGC/ASD/Grove_2019/iPSYCHPGC_ASD_Nov2017"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_ASD_Grove_2019_GR_p005.rds"))

# BP

fn <- "/binder/common/public_data/PGC/BP_2018/BDvsCONT.sumstats"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_BPD_GR_p005.rds"))

# MDD

library(biomaRt)
ensembl <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37")

fn <- "/binder/common/public_data/PGC/MDD_2019/PGC_UKB_depression_genome-wide_ACGT.txt"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("MarkerName", "P"))
gwas <- gwas[P <= 0.05]

gwas.mdd.biomart <- getBM(attributes = c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
              filters = "snp_filter", values = gwas$MarkerName, mart = ensembl, uniqueRows = TRUE)

gwas <- left_join(gwas, gwas.mdd.biomart, by = c("MarkerName" = "SNP")) %>% 
  dplyr::select(SNP, CHR, POS, P)

colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenomicRanges::GRanges(seqnames = gwas$CHR,
                                  ranges = IRanges::IRanges(start = as.numeric(as.character(gwas$POS)),
                                                            end = as.numeric(as.character(gwas$POS))),
                                  snp_id = gwas$SNP,
                                  p_value = gwas$P)  

gwas.gr <- gwas.gr[!duplicated(gwas.gr)] 
GenomeInfoDb::seqlevelsStyle(gwas.gr) <- "UCSC"

saveRDS(gwas.gr, file =  paste0(gwas.cluster.out.dir, "pgc_MDD_GR_p005.rds"))

# IBD

fn <- "/binder/common/public_data/IBD/EUR.IBD.gwas_info03_filtered.assoc"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "IBD_GR_p005.rds"))

# PD

fn <- "/binder/common/public_data/PD/daner_panic_3GDES_col1_12"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "PD_GR_p005.rds"))

# BMI

fn <- "/binder/common/public_data/BMI/Locke_et_al/SNP_gwas_mc_merge_nogc.tbl.uniq"

gwas.locke <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "p"))

fn <- "/binder/common/public_data/BMI/ENGAGE/ENGAGE1000G_BMI.txt"

gwas <- fread(fn) #, header = T, stringsAsFactors = F, select = c("SNP", "chromosome", "position"))
gwas <- gwas %>% dplyr::select(SNP, CHR = chromosome, POS = position)

gwas.merged <- dplyr::left_join(gwas.locke, gwas)

colnames(gwas.merged) <- c("SNP", "P", "CHR", "POS")

gwas.gr <- GenerateGrangesObject(gwas.merged, 
                                 ofile = paste0(gwas.cluster.out.dir, "BMI_GR_p005.rds"))
