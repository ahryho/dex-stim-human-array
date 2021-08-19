# Flip SNPs with negative effect size in PRS-CS output

# Read arguments
args   <- commandArgs(trailingOnly=TRUE)
trait  <- args[1]
cohort <- args[2]
prs.output.dir <- args[3]

# Function

flipRisk <- function(risks) # flip negative effect sizes to positive values
  {
  n <- ncol(risks)
  risks$V4 <- as.character(risks$V4)
  risks$V5 <- as.character(risks$V5)

  risks["FLIP"] <- F
  if(any(risks$V6<0)) risks[which(risks$V6<0),]$FLIP <- T

  risks$TEMP <- risks$V4
  if(any(risks$FLIP))
    {
    risks[which(risks$FLIP),]$V4 <- risks[which(risks$FLIP),]$V5
    risks[which(risks$FLIP),]$V5 <- risks[which(risks$FLIP),]$TEMP
    risks[which(risks$FLIP),]$V6 <- -risks[which(risks$FLIP),]$V6
    }
  risks <- risks[,1:n]
  return(risks)
  }


# Filter and flip

effects <- read.table(paste0(prs.output.dir, "/", trait, "/", cohort, "_", trait, "_all_chromosomes.txt"))

effects <- flipRisk(effects)

write.table(effects,
    paste0(prs.output.dir, "/", trait, "/", cohort, "_", trait, "_all_chromosomes.flipRisk"),
    r = F, col.names = F, qu = F, sep = "\t")
