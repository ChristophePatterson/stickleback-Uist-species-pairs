# Output file location
library(tidyverse)
library(patchwork)
args <- commandArgs(trailingOnly = TRUE)

# set path
my_bins <- args[1]
# my_bins <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi"
print(my_bins)

# read in data
sliding_wd <- as_tibble(read.csv(paste0(my_bins, ".csv"), header = T))

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

chr <- chr[chr$GenBank.seq.accession%in%unique(sliding_wd$scaffold),]
chr$Sequence.name <- gsub("chr", "", chr$Sequence.name)
# chr <- chr[order(chr$Seq.length, decreasing = T),]
chr$bi.col <- rep(1:2, length.out = dim(chr)[1])

sliding_wd$chr <- chr$Sequence.name[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]
sliding_wd$bi.col <- chr$bi.col[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]

## Rename Fst column to code is usualble across results
colnames(sliding_wd)[grep("Fst",colnames(sliding_wd))] <- "Fst"
colnames(sliding_wd)[grep("dxy",colnames(sliding_wd))] <- "dxy"

## Calculate summary stats
fst_stats <- sliding_wd %>%
  summarise(mn_fst = mean(Fst_anad_resi, na.rm = T),
            md_fst = median(Fst_anad_resi, na.rm = T),
            sd_fst = sd(Fst_anad_resi, na.rm = T),
            mx_fst = max(Fst_anad_resi, na.rm = T),
            mn_anad_pi = mean(pi_anad, na.rm = T),
            mn_resi_pi = mean(pi_anad, na.rm = T))

write.table(x = fst_stats, file = paste0(my_bins, "_sum_stats.txt"))

## Within ecotype genetic diversity
pi.plot <- ggplot(sliding_wd) +
  geom_point(aes(pi_resi, pi_anad), alpha = 0.25) +
  theme_classic() +
  facet_wrap(~chr)

ggsave(paste0(my_bins, "_pi_anad_vs_pi_resi.png"), pi.plot, width = 12, height = 12)

## Caclualte rolling average
num_windows <- 8
Fst_value <- 0.2
sliding_wd$rolling_av <- stats::filter(sliding_wd$Fst, filter = rep(1/3, num_windows), sides = 2)
## Convert to binary
sliding_wd$rolling_av <- as.numeric(sliding_wd$rolling_av>=Fst_value)
sliding_wd$rolling_av[sliding_wd$rolling_av==0] <- NA
# Set chr to order
chr$Sequence.name
sliding_wd$chr[1:10]
sliding_wd$chr <- factor(sliding_wd$chr, levels = chr$Sequence.name)
sliding_wd$chr[1:10]
## p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
##   #geom_point(show.legend = F, size = 0.2) +
##   geom_line(aes(y = Fst, x = mid, , col = as.factor(MF), fill = as.factor(bi.col)), linewidth = 0.2) +
##   # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
##   scale_color_manual(values = c("grey60", "black")) +
##   scale_fill_manual(values = c("grey60", "black")) +
##   ylim(c(0, max(sliding_wd$Fst))) +
##   facet_grid(chr~., scale = "free_x", space = "free_x") +
##   theme_classic() +
##   theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
##   ggtitle(basename(my_bins))
## p
## ggsave(filename = paste0(my_bins, "_line.pdf"), p, width = 10, height = 30)
## ggsave(filename = paste0(my_bins, "_line.png"), p, width = 10, height = 30)


p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  ylim(c(0, max(sliding_wd$Fst))) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(basename(my_bins))
p
ggsave(filename = paste0(my_bins, ".pdf"), p, width = 10, height = 30)
ggsave(filename = paste0(my_bins, ".png"), p, width = 10, height = 30)

p <- ggplot(sliding_wd, aes(x = mid, y = Fst, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  ylim(c(0, max(sliding_wd$Fst))) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(basename(my_bins))
p
ggsave(filename = paste0(my_bins, "_horizontal.pdf"), p, width = 30, height = 10)
ggsave(filename = paste0(my_bins, "_horizontal.png"), p, width = 30, height = 10)

ggsave(filename = paste0(my_bins, "_horizontal_mini.pdf"), p, width = 20, height = 4)
ggsave(filename = paste0(my_bins, "_horizontal_mini.png"), p, width = 20, height = 4)

## Plot specic regions

regions <- data.frame(chr = c("I"), start = c(26000000), end = c(27500000))

q <- ggplot(sliding_wd[sliding_wd$chr==regions$chr[1]&sliding_wd$start>regions$start[1]&sliding_wd$start<regions$end[1],]) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  ylim(c(0, max(sliding_wd$Fst))) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(basename(my_bins)) 

ggsave(filename = paste0(my_bins, "-", regions$chr[1], "-", regions$start[1], "-", regions$end[1], ".pdf"), q, width = 10, height = 5)
ggsave(filename = paste0(my_bins, "-", regions$chr[1], "-", regions$start[1], "-", regions$end[1], ".png"), q, width = 10, height = 5)


p <- ggplot(sliding_wd, aes(x = mid, y = dxy, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = dxy, ymin = 0, x = mid, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  ylim(c(0, max(sliding_wd$Fst))) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none") +
  ggtitle(basename(my_bins))
p
ggsave(filename = paste0(my_bins, "_dxy_horizontal.pdf"), p, width = 30, height = 10)
ggsave(filename = paste0(my_bins, "_dxy_horizontal.png"), p, width = 30, height = 10)


p <- ggplot(sliding_wd, aes(x = Fst, y = dxy)) +
  geom_point(alpha = 0.3) +
  # ylim(c(-0.1, max(sliding_wd$Fst))) +
  theme_classic() +
  ggtitle(basename(my_bins))
p
ggsave(filename = paste0(my_bins, "_dxy_vs_Fst_.pdf"), p, width = 20, height = 10)
ggsave(filename = paste0(my_bins, "_dxy_vs_Fst_.png"), p, width = 20, height = 10)
