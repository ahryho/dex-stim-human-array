setwd("~/bio/code/mpip/dex-stim-human-array/")
source("code/util.R")

### Load eQTM results
### 
out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"
eqtm.dex.nom.fn   <- paste0(out.dir.pre, "eqtms/eqtm_cis_result_dex_beta.csv") 
eqtm.dex.nom.df   <- fread(eqtm.dex.nom.fn, 
                           col.names = c("CpG_ID", "ENSG_ID", "beta_eqtm", "t-stat_eqtm", "p-value_eqtm", "fdr_eqtm"))

### Load meQTL opposite group
### 
dir.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/"

meqtl.opposite.fc.df <- fread(paste0(dir.pre, "meqtl_opposite_fc_gr_df.csv"))
opposite.fc.grp.cpg  <- meqtl.opposite.fc.df$CpG_ID %>% unique()
           
### Get overlap eQTM and meQTLs
###                                
venn.rslt     <- GetVennPlt(meqtl.df = meqtl.opposite.fc.df, eqtm.df = eqtm.dex.nom.df, cbPal.col = "#999999")
venn.eqtm.df  <- eqtm.dex.nom.df[CpG_ID %in% venn.rslt$cpgs, .(CpG_ID, ENSG_ID)] 
venn.meqtl.df <- meqtl.opposite.fc.df[CpG_ID %in% venn.rslt$cpgs][, .(meQTL_ID, CpG_ID, SNP)] %>% unique()

cpg.ids  <- venn.rslt$cpgs
ensg.ids <- venn.eqtm.df$ENSG_ID %>% unique()
snp.ids  <- venn.meqtl.df$SNP %>% unique()

### Save lists
###
write.table(cpg.ids,
           "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/parmegene_biomax/opposite_fc_intersect_meqtl_eqtm_cpgs.csv", row.names = F, col.names = F, quote = F)

write.table(ensg.ids,
            "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/parmegene_biomax/opposite_fc_intersect_meqtl_eqtm_ensgs.csv", row.names = F, col.names = F, quote = F)

write.table(snp.ids,
            "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/parmegene_biomax/opposite_fc_intersect_meqtl_eqtm_snps.csv", row.names = F, col.names = F, quote = F)


### Load data
### 
meth.veh <- LoadMethylBeta("veh")
meth.dex <- LoadMethylBeta("dex")

gex.veh <- LoadGEX("veh")
gex.dex <- LoadGEX("dex")

bio.veh <- fread( "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/bio_mtrx_methyl_gex_veh.csv")
bio.dex <- fread( "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/bio_mtrx_methyl_gex_dex.csv")

snps    <- LoadGenotype()

pheno   <- LoadPheno()
sample.veh.ids <- pheno[pheno$Dex == 0, c("DNA_ID", "RNA_ID", "DNAm_ID")]
sample.veh.ids <- sample.veh.ids[match(sample.veh.ids$DNA_ID, colnames(meth.veh)[-1]), ]

sample.dex.ids <- pheno[pheno$Dex == 1, c("DNA_ID", "RNA_ID", "DNAm_ID")]
sample.dex.ids <- sample.dex.ids[match(sample.dex.ids$DNA_ID, colnames(meth.dex)[-1]), ]


### Subset baseline
### 
meth.sub.veh <- meth.veh[CpG_ID %in% cpg.ids]
gex.sub.veh  <- gex.veh[ENSG_ID %in% ensg.ids]
snps.sub     <- snps[SNP %in% snp.ids]

### Save 
### 
out.dir.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/parmegene_biomax/"

fwrite(meth.sub.veh, 
       paste0(out.dir.pre, "veh_analysis/methyl_beta_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(gex.sub.veh, 
       paste0(out.dir.pre, "veh_analysis/gex_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(snps.sub, 
       paste0(out.dir.pre, "veh_analysis/snp_mtrx.csv"),
       quote = F, row.names = F, sep = ";")

### Subset dex
### 
meth.sub.dex <- meth.dex[CpG_ID %in% cpg.ids]
gex.sub.dex  <- gex.dex[ENSG_ID %in% ensg.ids]


### Veh + dex mtrx
### 
meth.sub.veh.dex <- cbind(meth.sub.veh, meth.sub.dex[, -1])
colnames(meth.sub.veh.dex) <- c("CpG_ID", sample.veh.ids$DNAm_ID, sample.dex.ids$DNAm_ID)

fwrite(meth.sub.veh.dex, 
       paste0(out.dir.pre, "veh_dex_analysis/methyl_beta_mtrx_veh_dex.csv"),
       quote = F, row.names = F, sep = ";")

### 
gex.sub.veh.dex <- cbind(gex.sub.veh, gex.sub.dex[, -1])
colnames(gex.sub.veh.dex) <- c("CpG_ID", sample.veh.ids$RNA_ID, sample.dex.ids$RNA_ID)

fwrite(gex.sub.veh.dex, 
       paste0(out.dir.pre, "veh_dex_analysis/gex_mtrx_veh_dex.csv"),
       quote = F, row.names = F, sep = ";")

###
bio.dex.veh           <- fread(paste0(out.dir.pre, "veh_analysis/bio_mtrx_gex_veh.csv"))
bio.sub.veh           <- bio.dex.veh[Feature %in% c("DNAm_SV1", "DNAm_SV2", "DNAm_SV3", "V1", "V2", "V3")]
colnames(bio.sub.veh) <- c("Feature", sample.veh.ids$DNAm_ID)
bio.sub.dex           <- bio.dex[Feature %in% c("DNAm_SV1", "DNAm_SV2", "DNAm_SV3", "V1", "V2", "V3")]
colnames(bio.sub.dex) <- c("Feature", sample.dex.ids$DNAm_ID)
bio.sub.veh.dex       <- cbind(bio.sub.veh, bio.sub.dex[, -1])

fwrite(bio.sub.veh.dex, 
       paste0(out.dir.pre, "veh_dex_analysis/bio_mtrx_veh_dex_svs_only.csv"),
       quote = F, row.names = F, sep = ";")

### mapping for CpG and ENSG
### 
mapping.full.df  <- fread("~/bio/code/mpip/dex-stim-human-array/data/mapping/mapping_cpg_gene_ensg_full.csv" )
mapping.bioxm.df <- na.omit(rbind(mapping.full.df[CpG_ID %in% cpg.ids], mapping.full.df[Ensemble_ID %in% ensg.ids]))[,.(CpG_ID, Ensemble_ID)]

fwrite(mapping.bioxm.df, 
       paste0(out.dir.pre, "mapping_parmegene_cpg_ensg.csv"),
       quote = F, row.names = F, sep = ";")
