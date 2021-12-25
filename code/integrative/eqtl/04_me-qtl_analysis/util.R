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
