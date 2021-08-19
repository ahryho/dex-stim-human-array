library(ggplot2)
library(tidyverse)
library(plyr)
library(data.table)

eigenval.fn         <- "/Users/anastasiia_hry/bio/workspace/dex/dex-genotype/pca/Dex_genoData_SNPs_LDpruned.eigenval"  # dex_geno_pca.eigenval"
eigenvec.fn         <- "/Users/anastasiia_hry/bio/workspace/dex/dex-genotype/pca/plink.mds" #pca_plink.eigenvec" #dex_geno_pca.eigenvec"
eigenvec.fn         <- "/Users/anastasiia_hry/bio/workspace/dex/dex-genotype/pca/Dex_genoData_SNPs_LDpruned.eigenvec"
pheno.snp.fn        <- "/Users/anastasiia_hry/bio/datasets/pheno/DexStim_LE_080216.txt"
pheno.methyl.fn     <- "/Users/anastasiia_hry/bio/datasets/pheno/pheno_full_for_kimono.csv"
pheno.methyl.pcs.fn <- "/Users/anastasiia_hry/bio/datasets/methylation/pheno_with_pcs.csv"


eigenval     <- read.table(eigenval.fn, sep = " ", header = F)
eigenvec     <- fread(eigenvec.fn)
pheno.snp    <- read.table(pheno.snp.fn, header = T)
pheno.methyl <- read.csv2(pheno.methyl.fn)



# 1. Take only those individuals that are in the methylation data

pheno.methyl <- pheno.methyl[!is.na(pheno.methyl$DNAm_ID),]
pheno.ids <- unique(pheno.methyl$DNA_ID)#Individual)
# eigenvec  <- eigenvec[eigenvec$V1 %in% pheno.ids, 2:ncol(eigenvec)]
eigenvec  <- eigenvec[, 2:ncol(eigenvec)]
names(eigenvec)[1] <- "Individual"
names(eigenvec)[3:ncol(eigenvec)] <- paste0("PC", 1:(ncol(eigenvec)-1))

ggplot(data = eigenvec, aes(PC1, PC2)) + geom_point()                       

# 2. Eigenvalues: % of variance

pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval) * 100)
names(pve)[2] <- "pve" 

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


# 3. K-means clustering
           
kmeans.mdl <- kmeans(eigenvec[, -c(1)], centers = 2, nstart = 100)            
clusters   <- kmeans.mdl$cluster
pca        <- cbind(eigenvec,  clusters)

# 4. PCA Individual map

ggplot(data = pca, aes(PC1, PC2, col = as.factor(clusters))) + 
  geom_point() +
  theme(legend.position = "none") + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
  

# 5. Corrplot

tmp.df <- data.frame(pheno.methyl) %>% dplyr::select(Dex, Sex, Status, Age, BMI_D1,
                                                     V1, V2, V3, V4, V5,
                                                     PC1, PC2,
                                                     DNAm_SV1, DNAm_SV2, DNAm_SV3)

tmp.df$DNAm_SV1 <- as.numeric(tmp.df$DNAm_SV1)
tmp.df$DNAm_SV2 <- as.numeric(tmp.df$DNAm_SV2)
tmp.df$DNAm_SV3 <- as.numeric(tmp.df$DNAm_SV3)

cor.mtrx <- cor(tmp.df)
corrplot(cor.mtrx, method = "color", type = "upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         sig.level = 0.01) 


# 6. Add PC1 and PC2 to pheno methyl data as covariates

pheno.methyl <- left_join(pheno.methyl, pca[, 1:3], by = "Individual")

write.table(pheno.methyl, file =pheno.methyl.pcs.fn, col.names = T, row.names = F, quote = F, sep = ";")

