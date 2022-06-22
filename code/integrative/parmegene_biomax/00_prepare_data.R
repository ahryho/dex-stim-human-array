setwd("~/bio/code/mpip/dex-stim-human-array/")

### Load eQTM results
### 
out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"
eqtm.dex.nom.fn   <- paste0(out.dir.pre, "eqtms/eqtm_cis_result_dex_beta.csv") 
eqtm.dex.nom.df   <- fread(eqtm.dex.nom.fn, 
                           col.names = c("CpG_ID", "ENSG_ID", "beta_eqtm", "t-stat_eqtm", "p-value_eqtm", "fdr_eqtm"))

### Load meQTL opposite group
### 
dir.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/"

meqtl.opposite.fc.df             <- fread(paste0(dir.pre, "meqtl_opposite_fc_gr_df.csv"))
opposite.fc.grp.cpg <- meqtl.opposite.fc.df$CpG_ID %>% unique()
           
### Get overlap eQTM and meQTLs
###                                
venn.rslt     <- GetVennPlt(meqtl.df = opposite.fc.grp.cpg, eqtm.df = eqtm.dex.nom.df, cbPal.col = "#999999")
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
