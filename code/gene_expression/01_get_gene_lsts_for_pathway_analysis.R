library(data.table)

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"
rslt.dir.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/gene_expression/02_gene_lists_for_pathway_analysis/"

# ----
# Load data

# Load eQTMs

eqtm.veh.nom.fn   <- paste0(out.dir.pre, "eqtms/eqtm_cis_result_veh_beta.csv") 
eqtm.dex.nom.fn   <- paste0(out.dir.pre, "eqtms/eqtm_cis_result_dex_beta.csv") 
eqtm.delta.nom.fn <- paste0(out.dir.pre, "eqtms/eqtm_cis_result_delta_beta.csv") 

col.names       <- c("CpG_ID", "ENSG_ID", "beta_eqtm", "t-stat_eqtm", "p-value_eqtm", "fdr_eqtm")

eqtm.veh.nom.df   <- fread(eqtm.veh.nom.fn, col.names = col.names)
eqtm.dex.nom.df   <- fread(eqtm.dex.nom.fn, col.names = col.names)
eqtm.delta.nom.df <- fread(eqtm.delta.nom.fn, col.names = col.names)

eqtm.all.nom.df   <- rbind(eqtm.veh.nom.df, eqtm.dex.nom.df, eqtm.delta.nom.df)

# Load meQTL groups

meqtl.parallel.fc.grp.veh.dex.gr <- readRDS(file =  paste0(out.dir.pre, "meqtls/primary_meqtl_results/meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_veh_dex_gr.rds"))
meqtl.parallel.fc.grp.delta.gr   <- readRDS(file =  paste0(out.dir.pre, "meqtls/primary_meqtl_results/meqtl_overlaps_meqtl_level_cpg_anno_parallel_fc_grp_delta_gr.rds"))
meqtl.opposite.fc.gr             <- readRDS(file =  paste0(out.dir.pre, "meqtls/primary_meqtl_results/meqtl_overlaps_meqtl_level_cpg_anno_opposite_fc_grp_gr.rds"))

parallel.fc.grp.veh.dex.cpg <- list(CpG_ID = names(meqtl.parallel.fc.grp.veh.dex.gr@anno))
parallel.fc.grp.delta.cpg   <- list(CpG_ID = names(meqtl.parallel.fc.grp.delta.gr@anno))
opposite.fc.grp.cpg         <- list(CpG_ID = names(meqtl.opposite.fc.gr@anno))

# Load GEX IDs
# 
gex.ids.fn <- "~/bio/code/mpip/dex-stim-human-array/data/gene_expression/mapping_ilmn_ensg_gene.csv"
gex.ids    <- fread(gex.mtrx.fn)

# ----
# Extract list of genes
# 

GetGenes <- function(meqtl.df, eqtm.df, out.fn){
  meqtl.cpgs <- meqtl.df$CpG_ID %>% unique()
  eqtm.cpgs  <- eqtm.df$CpG_ID %>% unique()
  
  intersect.cpgs <- intersect(meqtl.cpgs, eqtm.cpgs)
  ensg.ids <- eqtm.df[CpG_ID %in% intersect.cpgs, ENSG_ID] %>% unique()
  
  write.table(ensg.ids, 
             file = out.fn,
             row.names = F, col.names = F, quote = F)
}

# ----
# Background 
 
ensg.ids  <- gex.ids$Ensemble_ID %>% unique() # 11'959

out.fn <- paste0(rslt.dir.pre, "ensg_lst_bkgr.csv")
write.table(ensg.ids, 
            file = out.fn,
            row.names = F, col.names = F, quote = F)

# ----
# Target

# Parallel FC: veh & dex
#
out.fn <- paste0(rslt.dir.pre, "ensg_lst_parallel_fc_grp_veh_dex.csv")
GetGenes(meqtl.df = parallel.fc.grp.veh.dex.cpg, eqtm.df = eqtm.veh.nom.df, out.fn = out.fn) # 11'492

# Parallel FC: delta
#
out.fn <- paste0(rslt.dir.pre, "ensg_lst_parallel_fc_grp_delta.csv") # 2'277
GetGenes(meqtl.df = parallel.fc.grp.delta.cpg, eqtm.df = eqtm.dex.nom.df, out.fn = out.fn)

# Opposite FC
#
out.fn <- paste0(rslt.dir.pre, "ensg_lst_opposite_fc_grp.csv") #1'451
GetGenes(meqtl.df = opposite.fc.grp.cpg, eqtm.df = eqtm.dex.nom.df, out.fn = out.fn)
