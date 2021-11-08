---
  title: "eQTL Table and Figure Generation"
author:
  - name: Karolina Worf
affiliation: 
  - Institute of Computational Biology, Helmholtz Zentrum München, Neuherberg 85764, Germany
- TUM School of Life Sciences Weihenstephan, Technical University of Munich, Freising 85354, Germany
email: karolina.worf@helmholtz-muenchen.de
header-includes:
  - \usepackage{float}
- \usepackage{graphicx, multicol}
date: "`r doc_date()`"
abstract: In this R markdown file, we will provide the quality control and plot generation of the raw, normalized, and batch corrected data generated in the Microarray Preprocessing Pipeline.
output: 
  BiocStyle::pdf_document:
  toc: true
number_sections: true
editor_options: 
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Quality Control of the Affymetrix Microarrays}
%\VignetteEncoding{UTF-8}
%\VignetteEngine{knitr::rmarkdown}
---
  
  ```{r knitrOptions, include = F}
library(BiocStyle)
require(knitr)
options(stringsAsFactors = T)
opts_chunk$set(echo = TRUE,
               fig.pos = "H",
               dpi = 300,
               warning = FALSE, 
               message = FALSE, 
               cache = TRUE, 
               cache.lazy = FALSE)
```

<p>
  **R version**: `r R.version.string`
<br />
  **Bioconductor version**: `r BiocInstaller::biocVersion()`
</p>
  
  
  # Load Packages
  
  ```{r}
library(viridis)
library(UpSetR)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(reshape2)
library(VennDiagram)
library(oligo)
library(dplyr)
library(rvcheck)
library(ggplotify)
library(venneuler)
library(ggforce)
#library(beadarray)
library(tidyr)  # to separate IDs
```

# Load data

First, we will load all cis-eQTL datasets and filter them for only significant hits (FDR < 0.05).

```{r}
# cases and controls cis-eQTLs genes
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_CASE_CTL_matrixEQTL_cis_results_genes_probeset.txt", header = T, stringsAsFactors = F)
cis.genes.all <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05

# cases and controls cis-eQTLs transcripts
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_CASE_CTL_matrixEQTL_cis_results_transcripts.txt", header = T, stringsAsFactors = F)
cis.transcripts.all <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05

# cases and controls cis-eQTLs exons
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_CASE_CTL_matrixEQTL_cis_results_exons.txt", header = T, stringsAsFactors = F)
cis.exons.all <- data[data$FDR < 0.05,   ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05

# only cases cis-eQTLs genes
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/02_pmb_CASE_only_matrixEQTL_cis_results_genes_probeset.txt", header = T, stringsAsFactors = F)
cis.genes.cases <- data[data$FDR < 0.05,   ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05

# only cases cis-eQTLs transcripts
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/02_pmb_CASE_only_matrixEQTL_cis_results_transcripts.txt", header = T, stringsAsFactors = F)
cis.transcripts.cases <- data[data$FDR < 0.05,   ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05

# only cases cis-eQTLs exons
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/02_pmb_CASE_only_matrixEQTL_cis_results_exons.txt", header = T, stringsAsFactors = F)
cis.exons.cases <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
rm(data)
```

## Load SNP, gene, transcript and exon positions

As we need the SNP, gene, transcript and exon position for all matrixEQTL runs, we will load them in the beginning, so that we do not have to repeat this step over and over again.

Load SNP positions and save the necessary information:
  
  ```{r}
# Load snp positions
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", header = F, stringsAsFactors = FALSE) # snp chr pos
colnames(snpspos) <- c("chromosome", "SNP", "distance", "position", "SNP_allele1", "SNP_allele2")
snpspos$SNP_location <- paste(snpspos$chromosome, ":", snpspos$position, sep = "")
snpspos <- snpspos[, c("SNP_location", "SNP", "SNP_allele1", "SNP_allele2")]
```

Load gene positions and save the necessary information:
  
  ```{r}
genes <- import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf")
genes <- genes[genes$gene_type == "protein_coding"]  # get only protein-coding genes
genes <- genes[genes$type == "gene"]  # get only genes
genes <- as.data.frame(genes)
genes <- genes[, c("gene_id", "gene_name", "seqnames", "start", "end", "strand")]
genes$Gene_location <- paste(genes$seqnames, ":", genes$start, "..", genes$end, sep = "")
genespos <- genes[, c("gene_id", "gene_name", "Gene_location", "strand")]
colnames(genespos) <- c("Gene", "Gene_symbol", "Gene_location", "Gene_strand")
rm(genes)
```

Load transcript positions and save the necessary information:
  
  ```{r}
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
transpos <- transcriptSummaries@featureData@data
transpos$location <- gsub("-", "..", transpos$location)
rm(transcriptSummaries)
```

Load exon positions and save the necessary information:
  
  ```{r}
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
exonspos <- exonSummaries@featureData@data
exonspos$location <- gsub("-", "..", exonspos$location)
rm(exonSummaries)
```


## Significant Hits

Then, we will check the number of found eQTLs and unique gene names in all datasets.

|           | Genes ALL | Genes CASES only | Transcripts ALL | Transcripts CASES only | Exons ALL | Exons CASES only |
  |-----------|-----------|------------------|-----------------|------------------------|-----------|------------------|
  | cis-eQTLs | `r length(cis.genes.all$FDR)` | `r length(cis.genes.cases$FDR)` | `r length(cis.transcripts.all$FDR)` | `r length(cis.transcripts.cases$FDR)` | `r length(cis.exons.all$FDR)` | `r length(cis.exons.cases$FDR)` |
  | unique gene IDs | `r length(unique(cis.genes.all$gene))` | `r length(unique(cis.genes.cases$gene))` | `r length(unique(cis.transcripts.all$gene))` | `r length(unique(cis.transcripts.cases$gene))` | `r length(unique(cis.exons.all$gene))` | `r length(unique(cis.exons.cases$gene))` |
  | unique SNP IDs  | `r length(unique(cis.genes.all$SNP))` | `r length(unique(cis.genes.cases$SNP))` | `r length(unique(cis.transcripts.all$SNP))` | `r length(unique(cis.transcripts.cases$SNP))` | `r length(unique(cis.exons.all$SNP))` |`r length(unique(cis.exons.cases$SNP))` |
  
  
  ## Generate Result Tables
  
  For some genes, transcripts and exons the IDs are not uniquely identifiable, e.g. the ID is a combination of multiple IDs separated by a vertical line: ID1|ID2|ID3. We first separate these merged IDs into own rows and then add important information to the significant eQTL hits. In the end we save the data into result tables (.txt files). 

### Gene cis-eQTLs

#### ALL (cases + controls)

```{r}
# separate IDs into extra rows
cis.genes.all <- separate_rows(cis.genes.all, gene, sep = "\\|")
# add info to gene eqtl data (cases vs. controls) and save it into a file
eqtls.gene <- merge(cis.genes.all, snpspos, by = "SNP", all.x = T)
eqtls.gene <- merge(eqtls.gene, genespos, by.x = "gene", by.y = "Gene", all.x = T)
eqtls.gene <- eqtls.gene[complete.cases(eqtls.gene), ]
eqtls.gene <- eqtls.gene[, c("SNP", "gene", "Gene_symbol", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
eqtls.gene$eQTL_type <- "cis" 
colnames(eqtls.gene) <- c("SNP", "Gene", "Gene_symbol", "Beta", "Test_statistic", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location",  "Gene_strand", "eQTL_type")
eqtls.gene <- eqtls.gene[order(eqtls.gene$FDR),]
head(eqtls.gene)
write.table(eqtls.gene, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", row.names = F, quote = F, dec = ".", sep = " ")
```

#### CASES only

```{r}
# separate IDs into extra rows
cis.genes.cases <- separate_rows(cis.genes.cases, gene, sep = "\\|")
# add info to gene eqtl data (only cases) and save it into a file
eqtls.gene.cases <- merge(cis.genes.cases, snpspos, by = "SNP", all.x = T)
eqtls.gene.cases <- merge(eqtls.gene.cases, genespos, by.x = "gene", by.y = "Gene", all.x = T)
eqtls.gene.cases <- eqtls.gene.cases[complete.cases(eqtls.gene.cases), ]
eqtls.gene.cases <- eqtls.gene.cases[, c("SNP", "gene", "Gene_symbol", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
eqtls.gene.cases$eQTL_type <- "cis" 
colnames(eqtls.gene.cases) <- c("SNP", "Gene", "Gene_symbol", "Beta", "Test_statistic", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2","Gene_location",  "Gene_strand", "eQTL_type")
eqtls.gene.cases <- eqtls.gene.cases[order(eqtls.gene.cases$FDR),]
head(eqtls.gene.cases)
write.table(eqtls.gene.cases, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", row.names = F, quote = F, dec = ".", sep = " ")
```


### Transcript cis-eQTLs

#### ALL (cases + controls)

```{r}
# separate IDs into extra rows
cis.transcripts.all <- separate_rows(cis.transcripts.all, gene, sep = "\\|")
# add info to transcript eqtl data (cases vs. controls) and save it into a file
eqtls.transcript <- merge(cis.transcripts.all, snpspos, by = "SNP", all.x = T)
eqtls.transcript <- merge(eqtls.transcript, transpos, by.x = "gene", by.y = "transcript_id", all.x = T)
eqtls.transcript <- merge(eqtls.transcript, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.transcript <- eqtls.transcript[complete.cases(eqtls.transcript), ]
eqtls.transcript <- eqtls.transcript[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id",  "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.transcript) <- c("SNP", "Transcript", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Transcript_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.transcript$eQTL_type <- "cis" 
eqtls.transcript <- eqtls.transcript[order(eqtls.transcript$FDR),]
head(eqtls.transcript)
write.table(eqtls.transcript, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", row.names = F, quote = F, dec = ".", sep = " ")
```

#### CASES only

```{r}
# separate IDs into extra rows
cis.transcripts.cases <- separate_rows(cis.transcripts.cases, gene, sep = "\\|")
# add info to transcript eqtl data (only cases) and save it into a file
eqtls.transcript.cases <- merge(cis.transcripts.cases, snpspos, by = "SNP", all.x = T)
eqtls.transcript.cases <- merge(eqtls.transcript.cases, transpos, by.x = "gene", by.y = "transcript_id", all.x = T)
eqtls.transcript.cases <- merge(eqtls.transcript.cases, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.transcript.cases <- eqtls.transcript.cases[complete.cases(eqtls.transcript.cases), ]
eqtls.transcript.cases <- eqtls.transcript.cases[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id", "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.transcript.cases) <- c("SNP", "Transcript", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Transcript_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.transcript.cases$eQTL_type <- "cis" 
eqtls.transcript.cases <- eqtls.transcript.cases[order(eqtls.transcript.cases$FDR),]
head(eqtls.transcript.cases)
write.table(eqtls.transcript.cases, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", row.names = F, quote = F, dec = ".", sep = " ")
```


### Exon cis-eQTLs

#### ALL (cases + controls)

```{r}
# separate IDs into extra rows
cis.exons.all <- separate_rows(cis.exons.all, gene, sep = "\\|")
# add info to exon eqtl data (cases vs. controls) and save it into a file
eqtls.exon <- merge(cis.exons.all, snpspos, by = "SNP", all.x = T)
eqtls.exon <- merge(eqtls.exon, exonspos, by.x = "gene", by.y = "exon_id", all.x = T)
eqtls.exon <- merge(eqtls.exon, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.exon <- eqtls.exon[complete.cases(eqtls.exon), ]
eqtls.exon <- eqtls.exon[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id", "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.exon) <- c("SNP", "Exon", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.exon$eQTL_type <- "cis" 
eqtls.exon <- eqtls.exon[order(eqtls.exon$FDR),]
head(eqtls.exon)
write.table(eqtls.exon, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", row.names = F, quote = F, dec = ".", sep = " ")
```

#### CASES only

```{r}
# separate IDs into extra rows
cis.exons.cases <- separate_rows(cis.exons.cases, gene, sep = "\\|")
# add info to exon eqtl data (only cases) and save it into a file
eqtls.exon.cases <- merge(cis.exons.cases, snpspos, by = "SNP", all.x = T)
eqtls.exon.cases <- merge(eqtls.exon.cases, exonspos, by.x = "gene", by.y = "exon_id", all.x = T)
eqtls.exon.cases <- merge(eqtls.exon.cases, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.exon.cases <- eqtls.exon.cases[complete.cases(eqtls.exon.cases), ]
eqtls.exon.cases <- eqtls.exon.cases[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id", "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.exon.cases) <- c("SNP", "Exon", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.exon.cases$eQTL_type <- "cis" 
eqtls.exon.cases <- eqtls.exon.cases[order(eqtls.exon.cases$FDR),]
write.table(eqtls.exon.cases, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", row.names = F, quote = F, dec = ".", sep = " ")
```


## Upset plots

### SNP level

Upset plots (like VENN) to compare the datasets with each other:
  
  ```{r}
eqtls.snps <- list(genes_all = eqtls.gene$SNP, genes_cases_only = eqtls.gene.cases$SNP, transcripts_all = eqtls.transcript$SNP, transcripts_cases_only = eqtls.transcript.cases$SNP, exons_all = eqtls.exon$SNP, exons_cases_only = eqtls.exon.cases$SNP)
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/01_upset_plot_cis_eqtls_snps.pdf", width = 20, height = 7, onefile = F)
upset(fromList(eqtls.snps), 
      sets = c("exons_all", "transcripts_all", "genes_all", "exons_cases_only", "transcripts_cases_only", "genes_cases_only"), 
      nsets = 6, 
      nintersects = 50, 
      mainbar.y.label = "# SNP's", 
      text.scale = 2, 
      keep.order = T,
      sets.bar.color = c("#A6D854", "#4292C6", "#FDAE6B", "#A6D854", "#4292C6", "#FDAE6B"))
dev.off()
```

### Gene level

```{r}
eqtls.genes <- list(genes_all = eqtls.gene$Gene_symbol, genes_cases_only = eqtls.gene.cases$Gene_symbol, transcripts_all = eqtls.transcript$Gene_symbol, transcripts_cases_only = eqtls.transcript.cases$Gene_symbol, exons_all = eqtls.exon$Gene_symbol, exons_cases_only = eqtls.exon.cases$Gene_symbol)
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/01_upset_plot_cis_eqtls_gene_level.pdf", width = 20, height = 7, onefile = F)
upset(fromList(eqtls.genes), 
      sets = c("exons_all", "transcripts_all", "genes_all", "exons_cases_only", "transcripts_cases_only", "genes_cases_only"), 
      nsets = 6, 
      nintersects = 60, 
      mainbar.y.label = "# genes", 
      text.scale = 2, 
      keep.order = T,
      sets.bar.color = c("#A6D854", "#4292C6", "#FDAE6B", "#A6D854", "#4292C6", "#FDAE6B"))
dev.off()
```


## Combined boxplots of top exon eQTL

Find top exon hit without gene hit => nr. 1 top hit with SNP rs505460 and exon ENSE00001231633.1 of gene NCALD

```{r}
# unoverlapping genes
no.overlap <- setdiff(eqtls.exon$Gene_symbol, eqtls.gene$Gene_symbol)
length(no.overlap) # 3490
# overlapping genes
overlap <- intersect(eqtls.exon$Gene_symbol, eqtls.gene$Gene_symbol)
length(overlap) # 40% (n = 613) of 860 unique eQTL genes have an exon hit, while 3744 eQTL genes are only found with exons
# plot 2 Bsp overlap vs. non
exon.only <- eqtls.exon[eqtls.exon$Gene_symbol %in% no.overlap, ]
# get top hit
exon.only[which.min(exon.only$FDR), ]

# load necessary mapping and expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
huex <- fread("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/HuEx-1_0-st-v2/HuEx-1_0-st-v2.na36.hg19.probeset.csv", data.table = F, skip = 22)
huex <- huex[huex$start != "---", ]  # remove ids with no start and stop information
gene.code <- as.data.frame(import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf"))
# get gene expression
gene.id <- gene.code[gene.code$gene_name == "NCALD" & gene.code$type == "gene", "gene_id"]
gene.exprs <- t(exprs(geneSummaries[featureNames(geneSummaries) %in% gene.id, ])) # get expression of gene
# get transcript expression (only transcripts including the interesting exon)
transcript.ids <- gene.code[gene.code$gene_name == "NCALD" & gene.code$type == "transcript" & gene.code$transcript_name %in% c("NCALD-201", "NCALD-011", "NCALD-013", "NCALD-008"), "transcript_id"] 
transcript.exprs <- as.data.frame(204)
for (i in transcript.ids) {
  expr <- t(exprs(transcriptSummaries[grep(i, featureNames(transcriptSummaries)), ]))
  transcript.exprs <- cbind(transcript.exprs, expr)
}
transcript.exprs$`204` <- NULL
# get exon expression (only interesting exon/s)
exon.ids <- c("ENSE00001231633.1", "ENSE00002113070.1", "ENSE00002097865.1", "ENSE00002129643.1") # just take the interesting exon
exon.exprs <- as.data.frame(204)
for (i in exon.ids) {
  expr <- t(exprs(exonSummaries[grep(i, featureNames(exonSummaries)), ]))
  exon.exprs <- cbind(exon.exprs, expr)
}
exon.exprs$`204` <- NULL
# merge all three data sets
gt <- merge(gene.exprs, transcript.exprs, by = "row.names")
gte <- merge(gt, exon.exprs, by.x = "Row.names", by.y = "row.names")
head(gte)
save(gte, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/NCALD_rs505460_gene_transcript_exon_expression.RData") 

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --snp rs505460 --recode --tab --out NCALD_SNPs"))
snp <- read.delim("NCALD_SNPs.ped", head = T)
pheno <- pData(geneSummaries)[, c("IID","Dx")]

exp.p <- merge(gte, pheno, by.x = "Row.names", by.y = "row.names")
exp.p.snp <- merge(exp.p, snp, by.x = "IID", by.y = "X1.1")
save(exp.p.snp, file <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/NCALD_gene_transcript_exon_expression_snp.rda")

# Generate plot
plot.data <- melt(exp.p.snp, id.vars = "A.A", measure.vars = c("ENSG00000104490.13", grep("ENS", colnames(exp.p.snp), value = T)))
plot.data <- plot.data[plot.data$A.A != "0 0",]
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/03_boxplots_NCALD_rs505460_gene_transcript_exon_expression_snp.pdf")
ggplot(na.omit(plot.data), aes(A.A, value)) +
  geom_boxplot() + facet_wrap(~variable, scales = "free") +
  theme(strip.text = element_text(size = 4.1), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab("rs505460") + ylab("expression")
dev.off()
```


## Distance plots

Calculate the distances between the SNPs and genes/transcripts/exons in the cis-eQTLs.

```{r, results='hide'}
### GENE eQTL ALL ###
# 1. separate SNP chromosome from position
eqtls.gene$SNP_chr <- lapply(eqtls.gene$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.gene$SNP_pos <- lapply(eqtls.gene$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.gene$Gene_chr <- lapply(eqtls.gene$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.gene$Gene_start <- lapply(eqtls.gene$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.gene$Gene_end <- lapply(eqtls.gene$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.gene[, c("SNP_chr", "Gene_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.gene <- as.data.frame(eqtls.gene)

### GENE eQTL CASES only ###
eqtls.gene.cases$SNP_chr <- lapply(eqtls.gene.cases$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.gene.cases$SNP_pos <- lapply(eqtls.gene.cases$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.gene.cases$Gene_chr <- lapply(eqtls.gene.cases$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.gene.cases$Gene_start <- lapply(eqtls.gene.cases$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.gene.cases$Gene_end <- lapply(eqtls.gene.cases$Gene_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.gene.cases[, c("SNP_chr", "Gene_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.gene.cases <- as.data.frame(eqtls.gene.cases)

### Transcript eQTL ALL ###
eqtls.transcript$SNP_chr <- lapply(eqtls.transcript$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.transcript$SNP_pos <- lapply(eqtls.transcript$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.transcript$Transcript_chr <- lapply(eqtls.transcript$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.transcript$Transcript_start <- lapply(eqtls.transcript$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.transcript$Transcript_end <- lapply(eqtls.transcript$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.transcript[, c("SNP_chr", "Transcript_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.transcript <- as.data.frame(eqtls.transcript)

### Transcript eQTL CASES only ###
eqtls.transcript.cases$SNP_chr <- lapply(eqtls.transcript.cases$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.transcript.cases$SNP_pos <- lapply(eqtls.transcript.cases$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.transcript.cases$Transcript_chr <- lapply(eqtls.transcript.cases$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.transcript.cases$Transcript_start <- lapply(eqtls.transcript.cases$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.transcript.cases$Transcript_end <- lapply(eqtls.transcript.cases$Transcript_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.transcript.cases[, c("SNP_chr", "Transcript_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.transcript.cases <- as.data.frame(eqtls.transcript.cases)

### Exon eQTL ALL ###
eqtls.exon$SNP_chr <- lapply(eqtls.exon$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.exon$SNP_pos <- lapply(eqtls.exon$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.exon$Exon_chr <- lapply(eqtls.exon$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.exon$Exon_start <- lapply(eqtls.exon$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.exon$Exon_end <- lapply(eqtls.exon$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.exon[, c("SNP_chr", "Exon_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.exon <- as.data.frame(eqtls.exon)

### Exon eQTL CASES only ###
eqtls.exon.cases$SNP_chr <- lapply(eqtls.exon.cases$SNP_location, function(x) strsplit(x, ":")[[1]][1])
eqtls.exon.cases$SNP_pos <- lapply(eqtls.exon.cases$SNP_location, function(x) strsplit(x, ":")[[1]][2])
eqtls.exon.cases$Exon_chr <- lapply(eqtls.exon.cases$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][1], "chr")[[1]][2])
eqtls.exon.cases$Exon_start <- lapply(eqtls.exon.cases$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][1])
eqtls.exon.cases$Exon_end <- lapply(eqtls.exon.cases$Exon_location, function(x) strsplit(strsplit(x, ":")[[1]][2], "\\.\\.")[[1]][2])
# check visually if chromosomes of SNP and gene are equal
chromosomes <- as.data.frame(eqtls.exon.cases[, c("SNP_chr", "Exon_chr")])
chromosomes[!duplicated(chromosomes),] 
eqtls.exon.cases <- as.data.frame(eqtls.exon.cases)


# 2. get distance between gene and SNP
# Gene eQTLS ALL
dist.gene <- unlist(lapply(1:length(eqtls.gene$SNP), function(x) if(as.numeric(eqtls.gene[x, "SNP_pos"]) < as.numeric(eqtls.gene[x, "Gene_start"])){ as.numeric(eqtls.gene[x, "SNP_pos"]) - as.numeric(eqtls.gene[x, "Gene_start"])} else if(as.numeric(eqtls.gene[x, "SNP_pos"]) > as.numeric(eqtls.gene[x, "Gene_end"])){as.numeric(eqtls.gene[x, "SNP_pos"]) - as.numeric(eqtls.gene[x, "Gene_end"])} else {0}))

# Gene eQTLS CASES only
dist.gene.cases <- unlist(lapply(1:length(eqtls.gene.cases$SNP), function(x) if(as.numeric(eqtls.gene.cases[x, "SNP_pos"]) < as.numeric(eqtls.gene.cases[x, "Gene_start"])){ as.numeric(eqtls.gene.cases[x, "SNP_pos"]) - as.numeric(eqtls.gene.cases[x, "Gene_start"])} else if(as.numeric(eqtls.gene.cases[x, "SNP_pos"]) > as.numeric(eqtls.gene.cases[x, "Gene_end"])){as.numeric(eqtls.gene.cases[x, "SNP_pos"]) - as.numeric(eqtls.gene.cases[x, "Gene_end"])} else {0}))

# Transcript eQTLS ALL
dist.transcript <- unlist(lapply(1:length(eqtls.transcript$SNP), function(x) if(as.numeric(eqtls.transcript[x, "SNP_pos"]) < as.numeric(eqtls.transcript[x, "Transcript_start"])){ as.numeric(eqtls.transcript[x, "SNP_pos"]) - as.numeric(eqtls.transcript[x, "Transcript_start"])} else if(as.numeric(eqtls.transcript[x, "SNP_pos"]) > as.numeric(eqtls.transcript[x, "Transcript_end"])){as.numeric(eqtls.transcript[x, "SNP_pos"]) - as.numeric(eqtls.transcript[x, "Transcript_end"])} else {0}))

# Transcript eQTLS CASES only
dist.transcript.cases <- unlist(lapply(1:length(eqtls.transcript.cases$SNP), function(x) if(as.numeric(eqtls.transcript.cases[x, "SNP_pos"]) < as.numeric(eqtls.transcript.cases[x, "Transcript_start"])){ as.numeric(eqtls.transcript.cases[x, "SNP_pos"]) - as.numeric(eqtls.transcript.cases[x, "Transcript_start"])} else if(as.numeric(eqtls.transcript.cases[x, "SNP_pos"]) > as.numeric(eqtls.transcript.cases[x, "Transcript_end"])){as.numeric(eqtls.transcript.cases[x, "SNP_pos"]) - as.numeric(eqtls.transcript.cases[x, "Transcript_end"])} else {0}))

# Exon eQTLS ALL
dist.exon <- unlist(lapply(1:length(eqtls.exon$SNP), function(x) if(as.numeric(eqtls.exon[x, "SNP_pos"]) < as.numeric(eqtls.exon[x, "Exon_start"])){ as.numeric(eqtls.exon[x, "SNP_pos"]) - as.numeric(eqtls.exon[x, "Exon_start"])} else if(as.numeric(eqtls.exon[x, "SNP_pos"]) > as.numeric(eqtls.exon[x, "Exon_end"])){as.numeric(eqtls.exon[x, "SNP_pos"]) - as.numeric(eqtls.exon[x, "Exon_end"])} else {0}))

# Exon eQTLS CASES only
dist.exon.cases <- unlist(lapply(1:length(eqtls.exon.cases$SNP), function(x) if(as.numeric(eqtls.exon.cases[x, "SNP_pos"]) < as.numeric(eqtls.exon.cases[x, "Exon_start"])){ as.numeric(eqtls.exon.cases[x, "SNP_pos"]) - as.numeric(eqtls.exon.cases[x, "Exon_start"])} else if(as.numeric(eqtls.exon.cases[x, "SNP_pos"]) > as.numeric(eqtls.exon.cases[x, "Exon_end"])){as.numeric(eqtls.exon.cases[x, "SNP_pos"]) - as.numeric(eqtls.exon.cases[x, "Exon_end"])} else {0}))


# 3. generate distance plot
### GENE
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/04_distance_plot_gene_eqtls.pdf", height = 5)
plot(density(dist.gene[dist.gene != 0 & dist.gene > -1000000 & dist.gene < 1000000]), xlab = "Distance of gene eQTL-variant to TSS and TES", ylab = "Density", main = "", col = viridis(3)[1], ylim = c(0, 8e-06))
lines(density(dist.gene.cases[dist.gene.cases != 0 & dist.gene.cases > -1000000 & dist.gene.cases < 1000000]), col = viridis(3)[2])
legend("topright", legend = c("genes all", "genes cases only"), fill = c(viridis(3)[1],viridis(3)[2]))
dev.off()

### Transcript
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/04_distance_plot_transcript_eqtls.pdf", height = 5)
plot(density(dist.transcript[dist.transcript != 0 & dist.transcript > -1000000 & dist.transcript < 1000000]), xlab = "Distance of transcript eQTL-variant to TSS and TES", ylab = "Density", main = "", col = viridis(3)[1], ylim = c(0, 6e-06))
lines(density(dist.transcript.cases[dist.transcript.cases != 0 & dist.transcript.cases > -1000000 & dist.transcript.cases < 1000000]), col = viridis(3)[2])
legend("topright", legend = c("transcripts all", "transcripts cases only"), fill = c(viridis(3)[1],viridis(3)[2]))
dev.off()

### Exon
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/04_distance_plot_exon_eqtls.pdf", height = 5)
plot(density(dist.exon[dist.exon != 0 & dist.exon > -1000000 & dist.exon < 1000000]), xlab = "Distance of exon eQTL-variant to ESS and EES", ylab = "Density", main = "", col = viridis(3)[1], ylim = c(0, 7e-06))
lines(density(dist.exon.cases[dist.exon.cases != 0]), col = viridis(3)[2])
legend("topright", legend = c("exons all", "exons cases only"), fill = c(viridis(3)[1], virids(3)[2]))
dev.off()

# all combined
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/04_distance_plot_gene_transcript_exon_eqtls.pdf", height = 5)
plot(density(dist.gene[dist.gene != 0 & dist.gene > -1000000 & dist.gene < 1000000]), xlab = "Distance to TSS/ESS and TES/EES", ylab = "Density", main = "", col = viridis(6)[1], ylim = c(0, 8e-06))
lines(density(dist.gene.cases[dist.gene.cases != 0 & dist.gene.cases > -1000000 & dist.gene.cases < 1000000]), col = viridis(6)[2])
lines(density(dist.transcript[dist.transcript != 0 & dist.transcript > -1000000 & dist.transcript < 1000000]), col = viridis(6)[3])
lines(density(dist.transcript.cases[dist.transcript.cases != 0 & dist.transcript.cases > -1000000 & dist.transcript.cases < 1000000]), col = viridis(6)[4])
lines(density(dist.exon[dist.exon != 0 & dist.exon > -1000000 & dist.exon < 1000000]), col = viridis(6)[5])
lines(density(dist.exon.cases[dist.exon.cases != 0]), col = viridis(6)[6])
legend("topright", legend = c("genes all", "genes  cases only", "transcripts all", "transcripts cases only", "exons all", "exons cases only"), fill = viridis(6))
dev.off()
```

|                          | Only SCZ | Only CTL | ALL | 
  |--------------------------|----------|----------|-----------|
  | significant cis-eQTLs    | `r length(eqtls.scz$gene)` | `r length(eqtls.ctl$gene)` | `r length(eqtls.all$gene)` | 
  | eQTLs lying on gene      | `r length(dist.scz[dist.scz == 0])` | `r length(dist.ctl[dist.ctl == 0])` | `r length(dist.all[dist.all == 0])` | 
  | min distance all         | `r min(dist.scz)` | `r min(dist.ctl)` | `r min(dist.all)` | 
  | min distance not on gene | `r min(dist.scz[dist.scz != 0])` | `r min(dist.ctl[dist.ctl != 0])` | `r min(dist.all[dist.all != 0])` | 
  | max distance             | `r max(dist.scz)` | `r max(dist.ctl)` | `r max(dist.all)` | 
  | mean distance            | `r mean(dist.scz)` | `r mean(dist.ctl)` | `r mean(dist.all)` | 
  | median distance          | `r median(dist.scz)` | `r median(dist.ctl)` | `r median(dist.all)` |
  
  
  ```{r}
eqtls.exon <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = TRUE)
eqtls.exon <- eqtls.exon[eqtls.exon$Gene_symbol == "NCALD", ]
head(eqtls.exon)
snps <- GRanges(seqnames = eqtls.exon$SNP_location)
genes <- GRanges(seqnames = eqtls.exon$Gene_location)

```

------
  TODO
------
  
  ## Manhattan Plots
  
  1. The easy way using the qqman package:
  
  ```{r}
# for all
all <- merge(me.all$cis$eqtls, snpspos, by.x = "snps", by.y = "snp", all.x = T)
all <- all[complete.cases(all), ]
#png("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/manhattan_plot_gene_cis_eqtls_all.png")
manhattan(all, chr = "chr", bp = "pos", snp = "snps", p = "pvalue")
#dev.off()
# for SCZ only
cases <- merge(me.scz$cis$eqtls, snpspos, by.x = "snps", by.y = "snp", all.x = T)
cases <- cases[complete.cases(cases), ]
#png("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/manhattan_plot_gene_cis_eqtls_cases.png")
manhattan(cases, chr = "chr", bp = "pos", snp = "snps", p = "pvalue")
#dev.off()
# for CTL only
controls <- merge(me.ctl$cis$eqtls, snpspos, by.x = "snps", by.y = "snp", all.x = T)
controls <- controls[complete.cases(controls), ]
#bitmap("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/manhattan_plot_gene_cis_eqtls_controls.png", "png16m")
manhattan(controls, chr = "chr", bp = "pos", snp = "snps", p = "pvalue")
#dev.off()
```


# Boxplot of top eQTL with gene expression

Best hit is an eQTL with SNP ID rs151095999 and gene symbol ADGRG7. Now, we will plot the gene expression of this gene in the different alleles. 

```{r}
#plot the top 2
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --snps rs1133328 --recodeA --tab --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/eqtl-gex"))
snp <- read.delim("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/eqtl-gex.raw",head=T)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_core_gene_summaries_annotated_combat_batch_corrected.RData")

snp <- snp[, c("IID", "rs1133328_T")]
exp <- as.data.frame(cbind(t(exprs(geneSummaries[fData(geneSummaries)[fData(geneSummaries)$SYMBOL == "EFHB", "PROBEID"], ])), geneSummaries@phenoData@data$Dx, as.character(geneSummaries@phenoData@data$IID)))
colnames(exp) <- c("EFHB", "Dx", "IID")
eqtl <- merge(snp, exp, by = "IID")
# remove NA's
eqtl <- eqtl[complete.cases(eqtl), ]
dex0 <- eqtl[eqtl$Dx == 1, ]
dex1 <- eqtl[eqtl$Dx == 2 | eqtl$Dx == 3 | eqtl$Dx == 4, ]

a <- cbind(as.numeric(as.character(dex0$EFHB)), as.character(dex0$rs1133328_T), rep("Controls",16))
b <- cbind(as.numeric(as.character(dex1$EFHB)), as.character(dex1$rs1133328_T), rep("Cases",143))
colnames(a) <- c("EFHB", "rs1133328", "dx")
colnames(b) <- c("EFHB", "rs1133328", "dx")
c <- as.data.frame(rbind(a,b))
c$EFHB <- as.numeric(as.character(c$EFHB))

pdf("ILMN_1805449_fe.pdf")
ggplot(na.omit(c), aes(rs1133328, EFHB, fill = factor(dx))) +
  geom_boxplot() + facet_wrap(~ dx) + theme_classic()
dev.off()
```

### Scatterplots

#### eQTLs vs. own GWAS

```{r}

# Load eQTLs
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results_old/01_pmb_SCZ_CTL_matrixEQTL.RData")
me.all.gene <- me.all
eqtls.gene <- me.all.gene$cis$eqtls[me.all.gene$cis$eqtls$FDR < 0.05, c("snps", "pvalue")] # filter for significant eQTLs
rm(me.all, me.all.gene)
# Load GWAS
gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc.logistic", header = T, stringsAsFactors = F, data.table = F, select =  c("SNP", "P"))
gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values
gwas.eqtl.snps <- merge(eqtls.gene, gwasResults, by.x = "snps", by.y = "SNP") # merge the two datasets (not overlapping SNPs will be removed from the dataset)
gwas.eqtl.snps$log_pvalue <- -log10(gwas.eqtl.snps$pvalue)
gwas.eqtl.snps$log_P <- -log10(gwas.eqtl.snps$P)
ggplot(gwas.eqtl.snps, aes(x = log_pvalue, y = log_P)) + 
  geom_point(alpha = 0.3) +
  ylab("-log10(P) of GWAS") +
  xlab("-log10(P) of eQTLs")
```

#### eqTLs vs. SCZ GWAS

```{r}
# Load GWAS
gwasResults <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics//SCZ2/scz2.snp.results.txt", header = T, stringsAsFactors = F, data.table = F, select = c("snpid", "p"))
gwas.eqtl.snps <- merge(eqtls.gene, gwasResults, by.x = "snps", by.y = "snpid") # merge the two datasets (not overlapping SNPs will be removed from the dataset)
gwas.eqtl.snps$log_pvalue <- -log10(gwas.eqtl.snps$pvalue)
gwas.eqtl.snps$log_p <- -log10(gwas.eqtl.snps$p)
ggplot(gwas.eqtl.snps, aes(x = log_pvalue, y = log_p)) + 
  geom_point(alpha = 0.3) +
  ylab("-log10(P) of GWAS") +
  xlab("-log10(P) of eQTLs")
```


# Session Info

```{r, eval = F}
sessionInfo()
```

