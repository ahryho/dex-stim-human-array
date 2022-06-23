library(data.table)
out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/"

meqtl.delta.fn           <- paste0(out.dir.pre, "me-qtl_cis_indp_rw_delta_fdr_005.csv")
ind.meqtl.delta.df       <- fread(meqtl.delta.fn)

opposite.fc.delta.df     <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/meqtl_opposite_fc_gr_delta_df.csv")
parallel.fc.delta.df   <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/meqtl_parallel_fc_gr_delta_df.csv")


delta.cpgs <- ind.meqtl.delta.df$CpG_ID %>% unique() # 2'769

delta.meqtls <- ind.meqtl.delta.df[, .(CpG_ID, SNP)]  %>% unique() # 3'147
delta.meqtls[["meQTL_ID"]] <- paste(delta.meqtls$SNP, delta.meqtls$CpG_ID, sep = "-")

opposite.fc.delta.rw.df <- opposite.fc.delta.df[meQTL_ID %in% delta.meqtls$meQTL_ID]
parallel.fc.delta.rw.df <- parallel.fc.delta.df[meQTL_ID %in% delta.meqtls$meQTL_ID]

table(delta.meqtls$meQTL_ID %in% unique(parallel.fc.delta.rw.df$meQTL_ID)) # 2'900
table(delta.meqtls$meQTL_ID %in% unique(opposite.fc.delta.rw.df$meQTL_ID)) # 247

unique(opposite.fc.delta.rw.df$CpG_ID) %>% length() #245
unique(opposite.fc.delta.rw.df$SNP) %>% length() # 247

unique(parallel.fc.delta.rw.df$CpG_ID) %>% length() # 2529
unique(parallel.fc.delta.rw.df$SNP) %>% length() # 2771

table(delta.meqtls$meQTL_ID %in% c(parallel.fc.delta.rw.df$meQTL_ID, opposite.fc.delta.rw.df$meQTL_ID))
table(delta.meqtls$meQTL_ID %in% c(parallel.fc.delta.df$meQTL_ID, opposite.fc.delta.df$meQTL_ID))

# Save

fwrite(opposite.fc.delta.rw.df,
       "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/indep_meqtl_opposite_fc_gr_delta_df.csv",
       quote = F, row.names = F, sep = ";")

fwrite(parallel.fc.delta.rw.df,
       "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/meqtl_parallel_and_opposite_fc_groups/indep_meqtl_parallel_fc_gr_delta_df.csv",
       quote = F, row.names = F, sep = ";")
