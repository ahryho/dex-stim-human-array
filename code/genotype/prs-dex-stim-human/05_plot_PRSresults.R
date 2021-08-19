library(ggplot2)
library(data.table)

prs.res.dir <- "/binder/mgp/datasets/2020_DexStim_Array_Human/PRS/"


cohorts.list <- c("MARS_DexStim")
traits.list  <- list.dirs(path = prs.res.dir, full.names = F, recursive = F)
traits.list 

result.tbl <- data.frame(FID = character(),
                         IID = character(),
                         score_avg = numeric(),
                         cohort = character(),
                         trait = character())

# Go through all cohorts and traits
for (cohort in cohorts.list){
  for (trait in traits.list){
    print(paste(cohort, trait))
    prs.sumstat.tmp <- fread(paste0(prs.res.dir, trait, "/", cohort, "_", trait, "_individual_rs_plink_out_no_flip.sscore"), 
                             header = TRUE)
    prs.sumstat.tmp <- cbind(prs.sumstat.tmp[, c("#FID", "IID", "SCORE1_AVG")], rep(cohort, nrow(prs.sumstat.tmp)), rep(trait, nrow(prs.sumstat.tmp)))
    colnames(prs.sumstat.tmp) <- colnames(result.tbl)
    result.tbl <- rbind(result.tbl, prs.sumstat.tmp)
  }
}

# Plot results
pdf(file = paste0(prs.res.dir, "PRS_noFlip.pdf"))
ggplot(result.tbl, aes(x = trait, y = score_avg)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
