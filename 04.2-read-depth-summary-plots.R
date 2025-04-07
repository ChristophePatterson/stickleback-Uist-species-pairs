## Load packages
.libPaths(c("/gpfs01/home/mbzcp2/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))
library(ggplot2)
library(patchwork)
library(ggrepel)

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
bam_QC$std_coverage <- as.numeric(gsub("X", "", x =  bam_QC$std_coverage))
bam_QC$reads <- as.numeric(gsub(",", "", x =  bam_QC$reads))
bam_QC$mapped_reads <- as.numeric(gsub(",", "", x =  bam_QC$mapped_reads))
bam_QC$dupl_reads <- as.numeric(gsub(",", "", x =  bam_QC$dupl_read))

## Plot
p <- ggplot(bam_QC) +
  geom_point(aes(mn_coverage, mapped_reads, col = Ecotype))
q <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red")
s <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, Ave_map_qc, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, Ave_map_qc, col = Ecotype), height = 0) +
  geom_hline(yintercept = 10, col = "red")
r <- ggplot(bam_QC) +
  geom_boxplot(aes(Population, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Population, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
t <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, dupl_reads/mapped_reads, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, dupl_reads/mapped_reads, col = Ecotype), height = 0)

p1 <- (p /q / s / r / t)
# Save
ggsave(filename = "Mapping_coverage.png", p1, width = 10, height = 20)

###  Filter out files with low coverage
min_cov <- 5
min_QC <- 10
bam_QC_hiQ <- bam_QC[bam_QC$mn_coverage>=min_cov&bam_QC$Ave_map_qc>=min_QC,]
# Write out file for bam files to take forward
writeLines(con = "HiQ_bam_files.txt", paste(bam_QC_hiQ$sample, bam_QC_hiQ$bam_file))