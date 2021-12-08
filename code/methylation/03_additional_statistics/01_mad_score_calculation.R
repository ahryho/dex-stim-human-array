out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/"

mtrx.eqtl.rsl.dex <- readRDS(paste0(out.dir.pre, "me-qtl_matrx_dex.RDS"))
hist(mtrx.eqtl.rsl.dex$cis$hist.bins)


library(dplyr)
library(data.table)
library(ggplot2)

output.eqtm.pre <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"

pheno.fn           <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
pheno              <- fread(pheno.fn, na.strings = c('#N/A', "NA"), dec = ",") %>% setDT()
pheno              <- pheno[Include == 1]    

methyl.mtrx.veh <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_veh.csv"))
methyl.mtrx.dex <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_dex.csv"))
methyl.mtrx <- data.frame(cbind(methyl.mtrx.veh, methyl.mtrx.dex[, -1]))

# ALL

mad.all <- apply(methyl.mtrx[, -1], 1, function(x) mad(x, constant = 1))
mad.all.df <- data.frame(cbind(CpG_ID = methyl.mtrx$CpG_ID, MAD = mad.all)) 

summary(mad.all)

# VEH

mad.veh    <- apply(methyl.mtrx.veh[, -1], 1, function(x) mad(x, constant = 1))
mad.veh.df <- data.frame(cbind(CpG_ID = methyl.mtrx.veh$CpG_ID, MAD = mad.veh)) 
mad.veh.df$MAD <- as.numeric(as.character(mad.veh.df$MAD))

summary(mad.veh.df)

ggplot(mad.veh.df, aes(x = MAD)) +
  geom_density(alpha = 1, size = 1) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        #  panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8)) +
  xlim(0, 0.3) +
  labs(title = "DNAm MAD Score Density Plot, Baseline")

mad.veh.df <- mad.veh.df[order(mad.veh.df$MAD, decreasing = F), ]
mad.veh.df$MAD <- as.numeric(mad.veh.df$MAD)
mad.veh.df[, "perc"] <- cumsum(mad.veh.df$MAD) / sum(mad.veh.df$MAD) * 100

table(mad.veh.df$perc >= 95)
summary(mad.veh.df[mad.veh.df$perc >= 95, "MAD"])

mad.veh.df[, "treatment"] <- "veh"

# DEX

mad.dex <- apply(methyl.mtrx.dex[, -1], 1,  function(x) mad(x, constant = 1))
mad.dex.df <- data.frame(cbind(CpG_ID = methyl.mtrx.dex$CpG_ID, MAD = mad.dex)) 
mad.dex.df$MAD <- as.numeric(as.character(mad.dex.df$MAD))

summary(mad.dex.df)

ggplot(mad.dex.df, aes(x = MAD)) +
  geom_density(alpha = 1, size = 1) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        #  panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8)) +
  xlim(0, 0.3) +
  labs(title = "DNAm MAD Score Density Plot, Dex")

mad.dex.df <- mad.dex.df[order(mad.dex.df$MAD, decreasing = F), ]
mad.dex.df$MAD <- as.numeric(mad.dex.df$MAD)
mad.dex.df[, "perc"] <- cumsum(mad.dex.df$MAD) / sum(mad.dex.df$MAD) * 100

table(mad.dex.df$perc >= 90)

mad.dex.df[, "treatment"] <- "dex"

mad.df <- rbind(mad.veh.df, mad.dex.df)

fwrite(mad.df, 
       "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/mad_score_dnam_beta_values.csv",
       quote = F, row.names = F, sep = ";")

# Check
mad.df <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/methylation/mad_score_dnam_beta_values.csv")
mad.dex.df <- mad.df[treatment == "dex"]
mad.veh.df <- mad.df[treatment == "veh"]
