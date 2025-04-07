## Load packages
library(ggplot2)
library(patchwork)

## Get command line argument
args <- commandArgs(trailingOnly = TRUE)
bam_QC <- args[1]
out_dir <- args[2]
## Read in file
bam_QC <- read.table(bam_QC, header = T)
# Set working directory
setwd(out_dir)
## Read in sample data
sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-03-28.csv")

## Merge data with sample data
bam_QC <- merge(bam_QC, sample_data[,c(1,6:15)], by.x = "sample", by.y = "individual", all.x = T)

## Convert to numerics
bam_QC$mn_coverage <- as.numeric(gsub("X", "", x =  bam_QC$mn_coverage))
bam_QC$reads <- as.numeric(gsub(",", "", x =  bam_QC$reads))

## Plot
p <- ggplot(bam_QC) +
  geom_point(aes(mn_coverage, reads, col = Ecotype))
q <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype)) +
  geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype), height = 0)
# Save
ggsave(filename = "Mapping_coverage.png", p + q, width = 10, height = 5)