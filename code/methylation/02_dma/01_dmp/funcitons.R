library(ggplot2)

TuneFcandPPlot <- function(df, pval.list, fc.list, plot.title = " "){
  est.df           <- data.frame(matrix(NaN, nrow = length(pval.list), ncol = length(fc.list)))
  colnames(est.df) <- fc.list
  rownames(est.df) <- pval.list
  
  for (p.thr in pval.list)
    for (delta.fc in fc.list){
      dmp.sign.df  <- df[abs(FC) > delta.fc & FDR < p.thr, ]
      est.df[as.character(p.thr), as.character(delta.fc)] <- nrow(dmp.sign.df)
    }
  
  est.df <- rownames_to_column(est.df, var = "p-value")
  est.df <- est.df %>% pivot_longer(cols = colnames(est.df)[-1], names_to = "FC", values_to = "Nr_DMPs")
  
  est.df$`p-value` <- as.numeric(est.df$`p-value`)
  est.df$FC        <- est.df$FC
  
  g.plot <- 
    ggplot(est.df, aes(`p-value`, Nr_DMPs, color = FC, group = FC)) + 
      geom_point() + geom_line() +
      theme( panel.background = element_blank(),
           plot.title = element_text(size = 10),
           axis.title = element_text(size = 10),
           axis.text.x = element_text(angle = 0, hjust = 1), 
           legend.position = "bottom") +
      labs(title = plot.title, y = "Number of DMPs", x = "FDR BH")
  
  return(list(plot = g.plot, data = est.df))
}

VolcanoPlot <- function(df, pval, fc, model = "SV"){
  dmps.sign.df <- df[abs(FC) > fc & FDR < pval,]
  
  nr.hypermethyl <- nrow(dmps.sign.df[FC > fc & FDR < pval,])
  nr.hypomethyl  <- nrow(dmps.sign.df[FC < fc & FDR < pval,])
  
  df[, "threshold"] <- as.factor(abs(df$FC) > fc & df$FDR < pval)
  cols <- c("TRUE" = "green", "FALSE" = "grey")

  # grob <- grobTree(textGrob(nr.hypermethyl, x = fc/2 + fc,  y = -log10(pval) + 1, hjust=0,
  #                          gp = gpar(col = "red", fontsize = 10, fontface = "italic")))
  
  ggplot(df, aes(y = -log10(FDR), x = FC, color = threshold)) +
    geom_point(alpha = .5, size = 1.2) +
    scale_colour_manual(values = cols) +
    geom_vline(xintercept = fc, colour = "#990000", linetype = "dashed") + 
    geom_vline(xintercept = -fc, colour = "#990000", linetype = "dashed") +
    geom_hline(yintercept = -log10(pval), colour = "#990000", linetype = "dashed") +
    theme( panel.background = element_blank(),
           plot.title = element_text(size = 10),
           axis.title = element_text(size = 10),
           axis.text.x = element_text(angle = 0, hjust = 1), 
           legend.position = "none") +
    labs(title = paste0("Volcano plot: ", model, " model with FDR <= ", pval, " and abs(FC) >= ", fc),
         x = "Fold Change", y = "-log10 FDR") +
    annotate(geom = "text", x = 7*fc,  y = -log10(min(df$FDR))/2, label = nr.hypermethyl, color = "red", fontsize = 10, fontface = "bold") +
    annotate(geom = "text", x = - 7*fc,  y = -log10(min(df$FDR))/2, label = nr.hypomethyl, color = "red", fontsize = 10, fontface = "bold")
}