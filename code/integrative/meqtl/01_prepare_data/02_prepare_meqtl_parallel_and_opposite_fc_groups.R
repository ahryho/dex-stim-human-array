out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/"

meqtl.all.full.df <- fread(paste0(out.dir.pre, "meqtl_cis_result_veh_dex_delta_fdr_005.csv"))

### Read meQTLs with opposite effects

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/"

meqtls.opposite.fc.grp.with.delta <- 
  fread(paste0(out.dir.pre, "dex_veh_meqtl_opposite_beta_group_with_delta.csv"))

meqtls.opposite.fc.grp.no.delta <- 
  fread(paste0(out.dir.pre, "dex_veh_meqtl_opposite_beta_9K.csv"))

meqtl.opposite.fc <- c(meqtls.opposite.fc.grp.with.delta$meQTL_ID, meqtls.opposite.fc.grp.no.delta$meQTL_ID)

###
###
meqtls.meqtls <- list(delta = meqtl.all.full.df[treatment == "delta", meQTL_ID], 
                      dex = meqtl.all.full.df[treatment == "dex", meQTL_ID], 
                      veh = meqtl.all.full.df[treatment == "veh", meQTL_ID])

#####
#####
meqtl.grp.overlaps.all       <- intersect(intersect(meqtls.meqtls$veh, meqtls.meqtls$delta), meqtls.meqtls$dex)
meqtl.grp.overlaps.delta.veh <- setdiff(intersect(meqtls.meqtls$veh, meqtls.meqtls$delta), meqtls.meqtls$dex)
meqtl.grp.overlaps.delta.dex <- setdiff(intersect(meqtls.meqtls$dex, meqtls.meqtls$delta), meqtls.meqtls$veh)
meqtl.grp.only.delta         <- setdiff(setdiff(meqtls.meqtls$delta, meqtls.meqtls$dex), meqtls.meqtls$veh)

meqtl.grp.overlaps.dex.veh <- setdiff(intersect(meqtls.meqtls$veh, meqtls.meqtls$dex), meqtls.meqtls$delta)
meqtl.grp.only.dex         <- setdiff(setdiff(meqtls.meqtls$dex, meqtls.meqtls$veh), meqtls.meqtls$delta)
meqtl.grp.only.veh         <- setdiff(setdiff(meqtls.meqtls$veh, meqtls.meqtls$dex), meqtls.meqtls$delta)

meqtls.opposite.fc.grp.veh.dex <- setdiff(c(meqtl.grp.overlaps.dex.veh, meqtl.grp.only.dex, meqtl.grp.only.veh), meqtls.opposite.fc.grp.no.delta$meQTL_ID)
meqtls.opposite.fc.grp.delta   <- setdiff(c(meqtl.grp.overlaps.all, meqtl.grp.overlaps.delta.veh, meqtl.grp.overlaps.delta.dex,  meqtl.grp.only.delta), meqtls.opposite.fc.grp.with.delta$meQTL_ID)

#####
#####
meqtl.overlap.all.df       <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.overlaps.all]
meqtl.overlap.delta.veh.df <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.overlaps.delta.veh]
meqtl.overlap.delta.dex.df <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.overlaps.delta.dex]
meqtl.only.delta.df        <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.only.delta]
meqtl.overlap.dex.veh.df   <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.overlaps.dex.veh]
meqtl.only.dex.df          <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.only.dex]
meqtl.only.veh.df          <- meqtl.all.full.df[meQTL_ID %in% meqtl.grp.only.veh]

meqtl.opposite.fc.grp.with.delta.df <- meqtl.all.full.df[meQTL_ID %in% meqtls.opposite.fc.grp.with.delta$meQTL_ID]
meqtl.opposite.fc.grp.no.delta.df   <- meqtl.all.full.df[meQTL_ID %in% meqtls.opposite.fc.grp.no.delta$meQTL_ID]

### Group 1
meqtl.parallel.fc.grp.veh.dex.df  <- rbind(meqtl.overlap.dex.veh.df, 
                                           meqtl.only.dex.df, 
                                           meqtl.only.veh.df)[!(meQTL_ID %in% meqtls.opposite.fc.grp.no.delta$meQTL_ID)] 

fwrite(meqtl.parallel.fc.grp.veh.dex.df,
       paste0(out.dir.pre, "meqtl_parallel_and_opposite_fc_groups/meqtl_parallel_fc_gr_veh_dex_df.csv"),
       quote = F, row.names = F, sep   = "\t")

### Group 2
meqtl.parallel.fc.grp.delta.df    <- rbind(meqtl.overlap.all.df, 
                                           meqtl.overlap.delta.veh.df, 
                                           meqtl.overlap.delta.dex.df, 
                                           meqtl.only.delta.df)[!(meQTL_ID %in% meqtls.opposite.fc.grp.with.delta$meQTL_ID)]

fwrite(meqtl.parallel.fc.grp.delta.df,
       paste0(out.dir.pre, "meqtl_parallel_and_opposite_fc_groups/meqtl_parallel_fc_gr_delta_df.csv"),
       quote = F, row.names = F, sep   = "\t")

### Group 3
meqtl.opposite.fc.df <- rbind(meqtl.opposite.fc.grp.with.delta.df, 
                              meqtl.opposite.fc.grp.no.delta.df)

fwrite(meqtl.opposite.fc.df,
       paste0(out.dir.pre, "meqtl_parallel_and_opposite_fc_groups/meqtl_opposite_fc_gr_df.csv"),
       quote = F, row.names = F, sep   = "\t")
