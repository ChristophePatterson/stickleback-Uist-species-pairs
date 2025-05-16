# Output file location
library(tidyverse)
library(patchwork)
args <- commandArgs(trailingOnly = TRUE)

# set path
my_bins <- args[1]
print(my_bins)

# read in data
sliding_wd <- as_tibble(read.csv(paste0(my_bins, ".csv"), header = T))

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T))

chr <- chr[chr$RefSeq.seq.accession%in%unique(sliding_wd$scaffold),]
# chr <- chr[order(chr$Seq.length, decreasing = T),]
chr$bi.col <- rep(1:2, length.out = dim(chr)[1])

sliding_wd$chr <- chr$Sequence.name[match(sliding_wd$scaffold, chr$RefSeq.seq.accession)]
sliding_wd$bi.col <- chr$bi.col[match(sliding_wd$scaffold, chr$RefSeq.seq.accession)]

## Rename Fst column to code is usualble across results
colnames(sliding_wd)[grep("Fst",colnames(sliding_wd))] <- "Fst"

## Caclualte rolling average
num_windows <- 8
Fst_value <- 0.2
sliding_wd$rolling_av <- stats::filter(sliding_wd$Fst, filter = rep(1/3, num_windows), sides = 2)
## Convert to binary
sliding_wd$rolling_av <- as.numeric(sliding_wd$rolling_av>=Fst_value)
sliding_wd$rolling_av[sliding_wd$rolling_av==0] <- NA
# Set chr to order
sliding_wd$chr <- factor(sliding_wd$chr, levels = chr$Sequence.name)


p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  # ylim(c(-0.1, max(sliding_wd$Fst))) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(my_bins)
p
ggsave(filename = paste0(my_bins, ".pdf"), p, width = 10, height = 30)
ggsave(filename = paste0(my_bins, ".png"), p, width = 10, height = 30)

p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  # ylim(c(-0.1, max(sliding_wd$Fst))) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(my_bins)
p
ggsave(filename = paste0(my_bins, "_horizontal.pdf"), p, width = 30, height = 10)
ggsave(filename = paste0(my_bins, "_horizontal.png"), p, width = 30, height = 10)

## Plot specic regions

regions <- data.frame(chr = c("chrI"), start = c(25500000), end = c(27000000))

q <- ggplot(sliding_wd[sliding_wd$chr==regions$chr[1]&sliding_wd$start>regions$start[1]&sliding_wd$start<regions$end[1],]) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  # ylim(c(-0.1, max(sliding_wd$Fst))) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(my_bins) 

ggsave(filename = paste0(my_bins, "-", regions$chr[1], "-", regions$start[1], "-", regions$end[1], ".pdf"), q, width = 10, height = 5)
ggsave(filename = paste0(my_bins, "-", regions$chr[1], "-", regions$start[1], "-", regions$end[1], ".png"), q, width = 10, height = 5)
