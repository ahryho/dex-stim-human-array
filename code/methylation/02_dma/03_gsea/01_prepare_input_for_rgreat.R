ToBed <- function(df, output.fn, is.save = T){
  bed <- df %>% dplyr::select(chrom = chr, chromStart = pos, chromEnd = pos, name = PROBE_ID)
  rownames(bed) <- bed$PROBE_ID
  
  if (is.save == T)
    write.table(bed,
                output.fn, 
                sep = " ", quote = F, row.names = F, col.names = F)
  bed
}