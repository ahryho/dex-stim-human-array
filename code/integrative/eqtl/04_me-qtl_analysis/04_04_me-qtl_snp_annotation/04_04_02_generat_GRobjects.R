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
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "NCBI"
  print(length(input_gr)) # print length
  # save GRanges object to an rds file
  saveRDS(input_gr, file = ofile)
  return(input_gr)
}


out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

# background snps (baseline + delta )
# Read MAF values for each SNP
snps.dir <- "~/bio/code/mpip/dex-stim-human-array/data/snps/imputed_qc/from_janine/qc/"
mafs <- fread(paste0(snps.dir, "dex_geno_maf_for_filtered_samples.afreq"), data.table = F)
colnames(mafs) <- c("CHR", "SNP", "REF", "ALT", "MAF", "OBS_CT")

# generate bins of MAF in 0.05 steps
boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
mafs$bin <- arules::discretize(mafs$ALT_FREQS, method = "fixed", breaks = boundaries, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# load background SNPs
snps <- fread(paste0(snps.dir, "Dex_genoData_SNPs.bim"), data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos

background.all.gr <- GenerateGrangesObjecteQTL(input = snps, 
                                           mdata = mafs, 
                                           ofile = paste0(out.dir.pre, "dex_geno_snps_with_maf_bins_gr.rds"))

ind.meqtl.veh.df   <- fread(meqtl.veh.fn, col.names = col.names)
ind.meqtl.dex.df   <- fread(meqtl.dex.fn, col.names = col.names)

snps.bkgr <- snps[snps$SNP %in% ind.meqtl.delta.df$SNP,]

# Generate background SNP data as GRanges object
background.gr <- GenerateGrangesObjecteQTL(input = snps.bkgr, 
                                           mdata = mafs, 
                                           ofile = paste0(out.dir.pre, "meqtl_dex_snps_with_maf_gr.rds"))
table(seqnames(background.gr))  # no chr X/Y/M (remove in dataset if included there)
length(background.gr$snp_id) 

# Generate merging file for eQTL data from background SNPs 
background <- as.data.frame(background.gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")


# ind.meqtl.delta.gr <- readRDS(paste0(out.dir.pre, "meqtls_snp_annotated_withChIPseeker_delta.rds"))
# ind.meqtl.delta.df <- as.data.frame(ind.meqtl.delta.gr@anno)[, c("seqnames", "start")] %>% mutate(seqnames = sub("chr", "", seqnames))
# ind.meqtl.delta.df["SNP"] <- rownames(ind.meqtl.delta.df)
# colnames(ind.meqtl.delta.df) <- c("CHR", "POS", "SNP")


meqtl.delta.fn <- paste0(out.dir.pre, "me-qtl_cis_indp_result_delta_fdr_005.csv")
ind.meqtl.delta.df <- fread(meqtl.delta.fn, col.names = col.names)

meqtl.delta.snp.gr <- GenerateGrangesObjecteQTL(input = ind.meqtl.delta.df,
                                       mdata = background,
                                       ofile = paste0(out.dir.pre, "meqtl_delta_snps_with_maf_gr.rds"))

# GWAS Summary Statistics

# gwas.cross.dis <- read.table(text = gsub(" ", "\t", readLines("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_10k_oct2019_v2.txt.daner.txt")), header = T)
gwas.cross.dis <- fread("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt", header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"))
colnames(gwas.cross.dis) <- c("CHR", "POS", "SNP", "P")
gwas.cross.dis.10k <- fread("~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_10k_oct2019_v2.txt.daner.txt", header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"))
colnames(gwas.cross.dis.10k) <- c("CHR", "POS", "SNP", "P")

gwas.gr <- GenerateGrangesObject(gwas.cross.dis, ofile = "~/bio/code/mpip/dex-stim-human-array/data/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds")

