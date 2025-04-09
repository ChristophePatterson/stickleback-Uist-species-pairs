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
# Retain samples from focal waterbodies
paired_sp_waterbodies <- c("DUIN", "OBSE", "LUIB", "CLAC")
bam_QC_hiQ_paired <- bam_QC[bam_QC$Waterbody%in%paired_sp_waterbodies,]

## Plot just samples of interest
p <- ggplot(bam_QC_hiQ_paired) +
  geom_point(aes(mn_coverage, mapped_reads, col = Ecotype))
q <- ggplot(bam_QC_hiQ_paired) +
  geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red")
s <- ggplot(bam_QC_hiQ_paired) +
  geom_boxplot(aes(Ecotype, Ave_map_qc, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, Ave_map_qc, col = Ecotype), height = 0) +
  geom_hline(yintercept = 10, col = "red")
r <- ggplot(bam_QC_hiQ_paired) +
  geom_boxplot(aes(Population, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Population, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
t <- ggplot(bam_QC_hiQ_paired) +
  geom_boxplot(aes(Ecotype, dupl_reads/mapped_reads, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, dupl_reads/mapped_reads, col = Ecotype), height = 0)

p2 <- (p /q / s / r / t)
# Save
ggsave(filename = "Mapping_coverage_paired_populations.png", p2, width = 10, height = 20)

outgroup_locations <- c("Portugal", "Nova_Scotia", "Quebec", "Maine", "Iceland")
top.cov.outgroup <- list()
for(loc in outgroup_locations){
  bam_tmp <- bam_QC_hiQ[bam_QC_hiQ$Region==loc,]
  top.cov.outgroup[[loc]] <- bam_tmp[(order(bam_tmp$mn_coverage, decreasing = T)[1:2]),]
  print(paste("Top coverage samples for ", loc, " are:", paste(top.cov.outgroup[[loc]]$sample), "with ", paste(top.cov.outgroup[[loc]]$mn_coverage)))
}

top.cov.outgroup <- do.call("rbind", top.cov.outgroup)
bam_QC_hiQ_paired <- rbind(bam_QC_hiQ_paired, top.cov.outgroup)
# bam_QC_hiQ_paired

# Write out file for bam files to take forward
write.table(bam_QC_hiQ_paired, file = "HiQ_bam_files.csv", 
            row.names = F, quote = F, col.names = TRUE, sep = ",")
writeLines(con = "HiQ_bam_files.txt", paste(bam_QC_hiQ_paired$sample, bam_QC_hiQ_paired$bam_file))