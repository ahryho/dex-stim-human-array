out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/"

meqtl.dex.fn   <- paste0(out.dir.pre, "me-qtl_cis_result_dex_fdr_005.csv")
meqtl.veh.fn   <- paste0(out.dir.pre, "me-qtl_cis_result_veh_fdr_005.csv")
meqtl.delta.fn <- paste0(out.dir.pre, "me-qtl_cis_result_delta_fdr_005.csv")

col.names <-  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")

meqtl.dex.df   <- fread(meqtl.dex.fn, col.names = col.names)
meqtl.veh.df   <- fread(meqtl.veh.fn, col.names = col.names)
meqtl.delta.df <- fread(meqtl.delta.fn, col.names = col.names)

### Delta

meqtl.cpg.snp.ass     <- paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_delta_ind_cpg_snp_associations.txt")
meqtl.cpg.snp.ass.lst <- fread(meqtl.cpg.snp.ass) %>% unique()
meqtl.delta.ind.rw.df <- left_join(meqtl.cpg.snp.ass.lst, meqtl.delta.df) %>% unique()

meqtl.delta.ind.rw.df[["fdr"]] <- p.adjust(meqtl.delta.ind.rw.df$`p-value`, method = "BH")

fwrite(meqtl.delta.ind.rw.df,
       paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_indp_rw_delta_fdr_005.csv"),
       quote = F, row.names = F, sep = ";")

### Baseline

meqtl.cpg.snp.ass     <- paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_veh_ind_cpg_snp_associations.txt")
meqtl.cpg.snp.ass.lst <- fread(meqtl.cpg.snp.ass, header = F, fill = T) %>% unique()
colnames(meqtl.cpg.snp.ass.lst) <- c("CpG_ID", "SNP")
meqtl.cpg.snp.ass.lst <- meqtl.cpg.snp.ass.lst[SNP != ""]
meqtl.veh.ind.rw.df <- left_join(meqtl.cpg.snp.ass.lst, meqtl.veh.df)

meqtl.veh.ind.rw.df[["fdr"]] <- p.adjust(meqtl.veh.ind.rw.df$`p-value`, method = "BH")

fwrite(meqtl.veh.ind.rw.df,
       paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_indp_rw_veh_fdr_005.csv"),
       quote = F, row.names = F, sep = ";")

### Dex

meqtl.cpg.snp.ass     <- paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_dex_ind_cpg_snp_associations.txt")
meqtl.cpg.snp.ass.lst <- fread(meqtl.cpg.snp.ass, header = F, col.names =  c("CpG_ID", "SNP"))  %>% unique()
meqtl.dex.ind.rw.df   <- left_join(meqtl.cpg.snp.ass.lst, meqtl.dex.df)

meqtl.dex.ind.rw.df[["fdr"]] <- p.adjust(meqtl.dex.ind.rw.df$`p-value`, method = "BH")

fwrite(meqtl.dex.ind.rw.df,
       paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_indp_rw_dex_fdr_005.csv"),
       quote = F, row.names = F, sep = ";")
