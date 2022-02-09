Sys.setenv(BIOMART_CACHE="~/bio/caches/biomaRt")
Sys.setenv(ANNOTATION_HUB_CACHE="~/bio/caches/AnnotationHub")

GetFullmeQTLdf <- function(meqtl.df, fdr.thr = 0.05){
  colnames(meqtl.df) <- c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")
  
  meqtl.fltr.df <- meqtl.df[fdr < fdr.thr, ] %>% unique()
  
  meqtl.full.df <- left_join(meqtl.fltr.df, snp.bim) %>% mutate(pos_snp = pos) %>% select(-c(pos, pos_morgans))
  meqtl.full.df <- left_join(meqtl.full.df, cpg.loc[, .(CpG_ID, chr, pos)]) %>% 
    mutate(pos_cpg = pos) %>% select(-pos)
  meqtl.full.df <- meqtl.full.df[, dist := pos_snp - pos_cpg]  
  meqtl.full.df <- meqtl.full.df[, strand := ifelse(dist < 0, "-", "+")]
  meqtl.full.df <- meqtl.full.df[, dist := abs(dist)]
  
  meqtl.full.df
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