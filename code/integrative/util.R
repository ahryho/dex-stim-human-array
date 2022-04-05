Sys.setenv(BIOMART_CACHE="~/bio/caches/biomaRt")
Sys.setenv(ANNOTATION_HUB_CACHE="~/bio/caches/AnnotationHub")

GetFullmeQTLdf <- function(meqtl.df, fdr.thr = 0.05){
  colnames(meqtl.df) <- c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")
  
  meqtl.fltr.df <- meqtl.df[fdr < fdr.thr, ] %>% unique()
  
  meqtl.full.df <- left_join(meqtl.fltr.df, snp.bim) %>% mutate(pos_snp = pos) %>% 
    dplyr::select(-c(pos, pos_morgans))
  meqtl.full.df <- left_join(meqtl.full.df, cpg.loc[, .(CpG_ID, chr, pos)]) %>% 
    mutate(pos_cpg = pos) %>% 
    dplyr::select(-pos)
  meqtl.full.df <- meqtl.full.df[, dist := pos_snp - pos_cpg]  
  meqtl.full.df <- meqtl.full.df[, strand := ifelse(dist < 0, "-", "+")]
  meqtl.full.df <- meqtl.full.df[, dist := abs(dist)]
  
  meqtl.full.df
}

GetFulleQTMdf <- function(df, cpg.loc = NULL, ensg.loc = NULL, fdr.thr = 0.05, treatment, is.dist  = F){
  colnames(df) <- c("CpG_ID", "ENSG_ID", "beta", "t-stat", "p-value", "fdr")
  
  df <- df[fdr < fdr.thr, ] %>% unique()
  
  # if (isTRUE(is.dist) | (cpg.loc != NULL & ensg.loc != NULL)){
  #   df <- left_join(df, ensg.loc)
  #   df <- left_join(df, cpg.loc[, .(CpG_ID, chr, pos_cpg)])
  # }
  
  df[["treatment"]] <- treatment
  df[["eQTM_ID"]]   <- paste(df$CpG_ID, df$ENSG_ID, sep = "-")
  
  # meqtl.full.df <- meqtl.full.df[, dist := pos_snp - pos_cpg]  
  # meqtl.full.df <- meqtl.full.df[, strand := ifelse(dist < 0, "-", "+")]
  # meqtl.full.df <- meqtl.full.df[, dist := abs(dist)]
  
  df[, .(eQTM_ID, CpG_ID, ENSG_ID, beta, `p-value`, fdr, treatment)] %>% setDT()
}

# function for Scatter Plot

GetScatterPlot2 <- function(df, selected.meqtl = NULL, fdr.thr = 0.05, plot.title = NULL){
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(df, aes(x = beta, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1.5, size = 1.2) +
    scale_x_continuous(labels = scientific) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) +
    labs( x = "Effect size (FC)", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

GetScatterPlot3 <- function(df, selected.meqtl = NULL, fdr.thr = 0.05, plot.title = NULL){
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(df, aes(x = methyl_change, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1.5, size = 1.2) +
    scale_x_continuous(labels = scientific) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) +
    labs( x = "Effect size (FC)", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

# Function for Boxplot

ToBed <- function(df, output.fn, is.save = T){
  bed <- df %>% dplyr::select(chrom = chr, chromStart = pos, chromEnd = pos, name = PROBE_ID)
  rownames(bed) <- bed$PROBE_ID
  
  if (is.save == T)
    write.table(bed,
                output.fn, 
                sep = " ", quote = F, row.names = F, col.names = F)
  bed
}

## Scatter plot and Boxplots for significant meQTLs -->
  
# func-getbetavalues
GetBetValuesDF <- function(methyl.beta.df, snp.df, selected.qtl, treatment){
  beta.values.df    <- methyl.beta.df[CpG_ID %in% selected.qtl$CpG_ID, -1]
  beta.values.df    <- rbind(beta.values.df, snp.df[SNP %in% selected.qtl$SNP, -1])
  beta.values.df    <- data.frame(t(beta.values.df)) %>% setDT()
  colnames(beta.values.df) <- c("CpG", "SNP")
  beta.values.df$SNP <- as.factor(beta.values.df$SNP)
  beta.values.df[SNP == 0, SNP := "AA"]; beta.values.df[SNP == 1, SNP := "AB"]; beta.values.df[SNP == 2, SNP := "BB"] 
  
  beta.values.df$treatment <- treatment
  beta.values.df
}

# Scatter Plot
GetScatterPlot <- function(meqtl.all.full.df, selected.meqtl, fdr.thr = 0.05, plot.title = NULL){
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(meqtl.all.full.df, aes(x = beta, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1.5, size = 1.2) +
    geom_label_repel(data = selected.meqtl,
                     aes(x = beta,
                         y = -log10(fdr),
                         label = meQTL_ID),
                     fontface = 'bold',
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.5, "lines"),
                     segment.color = 'grey50',
                     nudge_x = 0.05, nudge_y = 10, 
                     size = 3) +
    scale_x_continuous(labels = scientific) +
    # scale_y_continuous(trans = trans_reverser('log10')) +
    # labs(title = " ", y = "", x = " ") + 
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          #  panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) +
    labs( x = "Effect size (FC), MatrixEQTL", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

# BoxPlot
GetBoxPlot <- function(beta.values.df, selected.meqtl, fdr.thr = 0.05, plot.labels = c("AA", "AB", "BB"), plot.title = NULL){
  
  cbPalette <- c(  "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  beta.values.df %>%
    ggplot(aes(y = CpG, x = SNP, fill = treatment)) +
    geom_boxplot(width = 0.2, color = "black") +
    # scale_fill_viridis(discrete = TRUE, alpha = 0.5) +
    scale_x_discrete(labels = plot.labels) +
    theme( panel.background = element_blank(),
           plot.title = element_text(size = 10),
           axis.title = element_text(size = 10),
           axis.text.x = element_text(angle = 0, hjust = 0.5), 
           legend.position = "bottom", 
           legend.title = element_blank()) +
    labs(y = paste0(selected.meqtl$CpG_ID, "\nDNAm beta value"), 
         x = selected.meqtl$SNP,
         title = plot.title) +
    scale_fill_manual(values = cbPalette)
}

# BoxPlot main
ProcessGetBoxPlot <- function(methyl.beta.veh.df, methyl.beta.dex.df, snp.df, selected.meqtl, 
                              fdr.thr = 0.05, plot.title = NULL){
  
  beta.values.dex.df <- GetBetValuesDF(methyl.beta.dex.df, snp.df, selected.meqtl, "dex")
  beta.values.veh.df <- GetBetValuesDF(methyl.beta.veh.df, snp.df, selected.meqtl, "veh")
  beta.values.df     <- rbind(beta.values.dex.df, beta.values.veh.df)
  
  snp.ind.cnt.df <- beta.values.dex.df %>% count(SNP) # beta.values.df[ , .(count = count(SNP)), by = .(SNP, treatment)]
  snp.ind.cnt.df[, label := paste0(SNP, " (n = ", n, ")")]
  
  GetBoxPlot(beta.values.df, selected.meqtl, snp.ind.cnt.df$label, fdr.thr = 0.05, plot.title)
}

GetManhattanPlot <- function(df, fdr.thr, ylims, plot.title){
  
  new.df <- df %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(CHR_LEN = max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(CHR_CUM_POS = cumsum(as.numeric(CHR_LEN)) - CHR_LEN) %>%
    dplyr::select(-CHR_LEN) %>%
    # Add this info to the initial dataset
    left_join(df, ., by = c("CHR" = "CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BP_CUM_POS = BP + CHR_CUM_POS)
  
  axis.df <- new.df %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BP_CUM_POS) + min(BP_CUM_POS)) / 2)
  
  cbPalette <- c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  ggplot(new.df, aes(x = BP_CUM_POS, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
    scale_color_manual(values = rep(cbPalette, 2)) + 
    scale_x_continuous(label = axis.df$CHR, breaks = axis.df$center, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    labs( x = "Chromosome",
          y = "-log10(P)",
          title = plot.title) +
    theme( panel.background = element_blank(),
           # plot.margin = margin(4, 2, 1, 1, unit = "mm"),
           plot.title = element_text(size = 10),
           axis.line.x = element_blank(),
           axis.line.y = element_blank(),
           panel.border = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           axis.title = element_text(size = 8),
           axis.text.x = element_text(angle = 0, hjust = 0.5), 
           legend.position = "none", 
           legend.title = element_blank()) 
}

# Take overlaps / non-overlaps

get_all_overlaps <- function(delta_df, dex_df, veh_df){
  intersect(intersect(veh_df, delta_df), dex_df)
}

get_veh_delta_overlaps <- function(delta_df, dex_df, veh_df){
  intersect(intersect(veh_df, delta_df), dex_df)
}

annotate_with_chipseeker <- function(df, chr, pos_start, pos_end, term = "SNP"){
  
  gr   <- GenomicRanges::GRanges(seqnames = paste0("chr", df[[chr]]),
                                 ranges = IRanges::IRanges(start = as.numeric(as.character(df[[pos_start]])),
                                                           end = as.numeric(as.character(df[[pos_end]]))))
  names(gr) <- df[[term]]
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  
  gr <- unique(gr)
  
  gr.anno <- annotatePeak(unique(gr), 
                          TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                          annoDb = "org.Hs.eg.db")
}
