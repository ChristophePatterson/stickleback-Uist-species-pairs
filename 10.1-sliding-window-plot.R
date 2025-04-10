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

colnames(sliding_wd)[grep("Fst",colnames(sliding_wd))] <- "Fst"

p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
  geom_line(show.legend = F, size = 0.2) +
  scale_color_manual(values = c("black", "grey60")) +
  ylim(c(0, max(sliding_wd$Fst))) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines')) +
  ggtitle(my_bins)

ggsave(filename = paste0(my_bins, ".pdf"), p, width = 30, height = 10)

