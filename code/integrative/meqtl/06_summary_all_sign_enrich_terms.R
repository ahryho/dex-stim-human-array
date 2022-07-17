library(data.table)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(RColorBrewer)

dir.pre     <- "~/bio/code/mpip/dex-stim-human-array/"
out.dir.pre <- paste0(dir.pre, "output/data/integrative/matrixEQTL/meqtls/")

ind.meqtl.delta.df      <- fread(paste0(out.dir.pre, "region_wise_independent_snps/me-qtl_cis_indp_rw_delta_fdr_005.csv"))
delta.meqtl.cpgs        <- ind.meqtl.delta.df$CpG_ID %>% unique() 

col.names <- c("CpG_ID", "Type", "SNP")

p.val.thrsh <- 0.05

# CpG enrichment

## Realtion to island

delta.meqtl.cpg.anno.df  <- fread(paste0(out.dir.pre, "region_wise_independent_snps/meqtls_cpg_annotated_delta.csv"), sep = "\t")
gen.loc.enrich.perm.rslt <- fread(paste0(dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_cpgs_relation_to_island_enrichment_perm.csv"))

types          <- gen.loc.enrich.perm.rslt[p_val_emp < p.val.thrsh, Feature]
cpg.isl.df     <- delta.meqtl.cpg.anno.df[Relation_to_Island %in% types, .(Name, Relation_to_Island)]
cpg.isl.df     <- left_join(cpg.isl.df, ind.meqtl.delta.df[, .(CpG_ID, SNP)], by = c("Name" = "CpG_ID"))
colnames(cpg.isl.df) <- col.names
# cpg.isl.df.lst <- split(cpg.isl.df, cpg.isl.df$Relation_to_Island)

## In relation to gene
delta.meqtl.snp.anno.rds  <- readRDS(paste0(out.dir.pre, "region_wise_independent_snps/meqtls_cpg_annotated_withChIPseeker_delta_gr.rds")) %>% 
  as.data.frame() %>% setDT()

gen.loc.enrich.perm.rslt <- fread(paste0(dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_cpgs_chipseeker_enrichment_perm_delta_vs_veh.csv"))

types           <- gen.loc.enrich.perm.rslt[p_val_emp < p.val.thrsh, Feature]
cpg.gene.loc.df <- delta.meqtl.snp.anno.rds[annotation %in% types, .(CpG_ID, annotation)]
cpg.gene.loc.df <- left_join(cpg.gene.loc.df, ind.meqtl.delta.df[, .(CpG_ID, SNP)])
colnames(cpg.gene.loc.df) <- col.names

# Overlap with GR-DMPs
# 
dmps.cpgs <- fread(paste0(dir.pre, "output/data/methylation/02_dmp/dmps_fdr01_fc02_anno_full.csv"), select = "PROBE_ID", col.names = "CpG_ID")

dmps.intersect <- intersect(ind.meqtl.delta.df$CpG_ID, dmps.cpgs$CpG_ID)
dmps.intersect.df <- data.frame(CpG_ID = dmps.intersect, Type = "GR-DMPs")
dmps.intersect.df <- left_join(dmps.intersect.df, ind.meqtl.delta.df[, .(CpG_ID, SNP)])

# ChromHMM enrichment

epigenome <- fread("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/epigenomes.tsv",
                   col.names = c("EID", "ORDER", "GROUP", "MNEMONIC", "NAME"))

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/"

delta.meqtl.snp.chromhmm.anno.df <- fread(paste0("~/bio/code/mpip/dex-stim-human-array/", "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtls_snp_chromhmm_annotated_delta.csv"), sep = "\t")

## Blood

tissue <- "blood"
enrichment.rslts.fn <- paste0(out.dir.pre, "meqtl_snps_chromHMM_", tissue, "/", list.files(paste0(out.dir.pre, "meqtl_snps_chromHMM_", tissue)))

chromhmm.blood.df   <- lapply(enrichment.rslts.fn, function(x){
  eid   =  sub(".csv", "", sub(".*_", "", x))
  ename = epigenome[EID == eid, NAME] 
  df    = read.csv2(x)
  return(data.frame(df, eid = eid, ename = ename))
}) %>% bind_rows() %>% setDT()

types    <- chromhmm.blood.df[p_val_emp < p.val.thrsh, .(state, eid, ename)]
blood.df <- delta.meqtl.snp.chromhmm.anno.df[annot.name %in% types$state][annot.code %in% types$eid]
blood.df <- left_join(blood.df, ind.meqtl.delta.df[, .(CpG_ID, SNP)])
blood.df[, Type := annot.name]#paste0(annot.name, "-", annot.code)]
blood.df <- blood.df[,.(CpG_ID, Type, SNP)] %>% unique()

# GWAS enrichment
# 


df <- rbind(cpg.isl.df, cpg.gene.loc.df, blood.df, dmps.intersect.df)
df.lst <- split(df[, c("Type", "SNP")], df$Type)
df.lst <- lapply(df.lst, function(df) df$SNP)

nb.cols <- length(unique(df$Type))
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

m <- as.data.frame(list_to_matrix(df.lst))

library(UpSetR)

upset(m, 
      sets = unique(df$Type), # names(df.lst),
      nsets = length(df.lst), 
      queries = list(list(query = intersects, params = list("Island", "Promoter (<=1kb)", "15_Quies"), color = "red"),
                     list(query = intersects, params = list("Island", "Promoter (<=1kb)", "15_Quies", "14_ReprPCWk"), color = "red"),
                     list(query = intersects, params = list("Island", "Promoter (<=1kb)", "15_Quies", "14_ReprPCWk", "5_TxWk"), color = "red"),
                     list(query = intersects, params = list("Island", "Promoter (<=1kb)", "15_Quies", "4_Tx", "5_TxWk"), color = "red")),
      nintersects = 73, 
      sets.bar.color = mycolors,
      mainbar.y.label = "Nr of intersections\nmeSNPs", 
      text.scale = 1.3, 
      keep.order = T,
      line.size = 0.1,
      mb.ratio = c(0.3, 0.7),
      # sets.bar.color = c( "#0072B2", "#009E73", "#E69F00"), 
      order.by = c("freq"), # "degree"), 
      group.by = "degree")


df.lst.interest <- list(`7_Enh` = df.lst$`7_Enh`,
                        `GR-DMPs` = df.lst$`GR-DMPs`,
                        `Promoter (<=1kb)` = df.lst$`Promoter (<=1kb)`,
                        Island = df.lst$Island)

upset(fromList(df.lst.interest), 
      sets = names(df.lst.interest),
      nsets = length(df.lst.interest), 
      nintersects = 50, 
      sets.bar.color = mycolors[1:4],
      mainbar.y.label = "Number of intersections, meSNPs", 
      text.scale = 1, 
      keep.order = T,
      line.size = 0.1,
      # sets.bar.color = c( "#0072B2", "#009E73", "#E69F00"), 
      order.by = "freq")
