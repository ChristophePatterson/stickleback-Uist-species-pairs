library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly=T)

plot.dir <- args[1]
# plot.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/pca/Anad_resi_fw/wndsize100000_wndslid50000/"
pca_mds_file <- gsub(".txt", "", args[2])
# pca_mds_file <- gsub(".txt", "", "sliding-window_pca_wndsize100000_wndslid50000.txt")

pca.comp.df <- read_csv(paste0(plot.dir, pca_mds_file, ".txt"))
# pca.comp.df <- read_csv("/gpfs01/home/mbzcp2/data/sticklebacks/results/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/pca/Anad_resi_fw/wndsize100000_wndslid50000/sliding-window_pca_wndsize25000_wndslid5000.txt")

# Read in chromosome details
chr <- chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T))
# Remove smaller scafs and mito
chr <- chr[!chr$Chromosome.name%in%c("Un","MT"),]
# Order Chromsome name to be in order 
chr$Chromosome.name <- factor(chr$Chromosome.name, levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

# Order chr in order in pca
pca.comp.df$chr <- chr$Chromosome.name[match(pca.comp.df$chr, chr$RefSeq.seq.accession)]
pca.comp.df$chr <- factor(pca.comp.df$chr, levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

# Create column each window for eaiser subset
pca.comp.df$windowname <- paste(pca.comp.df$chr, pca.comp.df$start, pca.comp.df$end, sep = "-")

## Add in sample data
sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
pca.comp.df$Population <- sample_data$Population[match(pca.comp.df$sample, sample_data$individual)]
pca.comp.df$Ecotype <- sample_data$Ecotype[match(pca.comp.df$sample, sample_data$individual)]
pca.comp.df$Waterbody <- sample_data$Waterbody[match(pca.comp.df$sample, sample_data$individual)]

# Remove stream
pca.comp.df$Ecotype[pca.comp.df$Ecotype=="st"] <- "fw"
## Change level
pca.comp.df$Ecotype <- factor(pca.comp.df$Ecotype, levels = c("anad", "resi", "fw")) 

# Transform PCA so that the axis is also segregating populations in the same direction across all windows
## Copy over PCA and MDS data to new scaled columns
pca.comp.df$MDS1_scaled <- pca.comp.df$MDS1
pca.comp.df$MDS2_scaled <- pca.comp.df$MDS2
pca.comp.df$PCA1_scaled <- pca.comp.df$PCA1
pca.comp.df$PCA2_scaled <- pca.comp.df$PCA2

## Convert to numeric
pca.comp.df$end <- as.numeric(pca.comp.df$end)

# Define the columns you want to check and potentially invert
scale_cols <- c("MDS1_scaled", "MDS2_scaled", "PCA1_scaled", "PCA2_scaled")

# Compute the sign for each window and ecotype == "anad"
signs <- pca.comp.df %>%
  filter(Ecotype == "anad") %>%
  group_by(windowname) %>%
  summarise(across(all_of(scale_cols), ~ sign(median(.x, na.rm = TRUE)), .names = "sign_{.col}"), .groups = "drop")

# Join sign info back to original data
pca.comp.df <- pca.comp.df %>%
  left_join(signs, by = "windowname") %>%
  mutate(across(all_of(scale_cols),
                ~ ifelse(get(paste0("sign_", cur_column())) == -1, -.x, .x))) %>%
  select(-starts_with("sign_"))  # remove helper columns

# Calculate max and min MDS
max_min <- pca.comp.df %>%
  group_by(windowname) %>%
  summarise(
    max.mds = abs(max(MDS1_scaled, na.rm = TRUE)),
    min.mds = abs(min(MDS1_scaled, na.rm = TRUE)),
    .groups = 'drop'
  )

## Calculate position of individual samples within max and min MDS
pca.comp.df <- pca.comp.df %>%
  left_join(max_min, by = "windowname") %>%
  mutate(
    MDS1_ratio = (MDS1_scaled + min.mds) / (max.mds + min.mds)
    )

## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))
# Save output
ggsave(paste0(plot.dir, pca_mds_file,"_pca12_point.png"), p/q, width = 40, height = 15)

## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))
# Save output
ggsave(paste0(plot.dir, pca_mds_file,"_pca12_line.png"), p/q, width = 40, height = 15)

# Break down by chromosome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)), name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file,"_Popsplit.png"), p, width = 40, height = 20)

## PLot MDS1 along genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds1.png"), p, width = 40, height = 20)

## PLot MDS1 along genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(chr~., scale = "free_x", space = "free_x", switch = "y") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds1_vert.png"), p, width = 30, height = 40)
ggsave(paste0(plot.dir, pca_mds_file, "_mds1_vert.pdf"), p, width = 30, height = 40)

## Plot MDS1 along genome and split by waterbody
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "mds_Popsplit.png"),p , width = 40, height = 20)

## Plot MDS1 along genome and split by waterbody
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "line_Popsplit.png"), p, width = 40, height = 20)
ggsave(paste0(plot.dir, pca_mds_file, "line_Popsplit.pdf"), p, width = 40, height = 20)

## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), MDS2_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))
# Save output
ggsave(paste0(plot.dir, pca_mds_file,"mds_line.png"), p/q, width = 40, height = 15)


### Zoomed in sections

regions <- data.frame(chr = c("I", "IX", "XI", "XXI"), start = c(25000000, 4500000, 5000000, 8000000), end = c(31000000, 10000000, 10000000, 15000000))

# Filer dataset to specific region
pca.comp.df.filt <- pca.comp.df %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

regions_plot <- ggplot(pca.comp.df.filt,
                       aes(as.numeric(end), MDS1_scaled, col = Population)) +
  #geom_segment(data = Venu_2022_Inv, aes(x = start, xend = end, y = min(pca.comp.df.filt$MDS1_scaled), yend = min(pca.comp.df.filt$MDS1_scaled)),
  #                                       col = "orange") +
  #geom_segment(data = Venu_2022_Inv_v5, aes(x = start, xend = end, y = max(pca.comp.df.filt$MDS1_scaled), yend = max(pca.comp.df.filt$MDS1_scaled)),
  #                                       col = "red") +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr,scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows.png"), regions_plot, width = 40, height = 15)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows.pdf"), regions_plot, width = 40, height = 15)

tile_plot <- ggplot(pca.comp.df.filt,
                       aes(as.numeric(end), sample, fill = MDS1_ratio, shape = Ecotype)) +
  geom_tile() +
  scale_fill_gradient2(low = "deepskyblue", mid = "orange" ,high = "darkgreen", midpoint=0.5) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5),
        panel.background = element_rect(fill = NA, color = "black"), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_specificWindows.png"), tile_plot, width = 40, height = 15)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_specificWindows.pdf"), tile_plot, width = 40, height = 15)

tile_plot_chrI <- ggplot(pca.comp.df[pca.comp.df$chr=="I"&pca.comp.df$start>=25900000&pca.comp.df$end<=26700000,],
                       aes(as.numeric(end), sample, fill = MDS1_ratio, shape = Ecotype)) +
  geom_tile() +
  scale_fill_gradient2(low = "deepskyblue", mid = "orange" ,high = "darkgreen", midpoint=0.5) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5),
        panel.background = element_rect(fill = NA, color = "black"), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_chI_inv.png"), tile_plot_chrI , width = 20, height = 15)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_chI_inv.pdf"), tile_plot_chrI , width = 20, height = 15)


## Inversion Identified by Venu 2022
Venu_2022_Inv <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Venu-et-al-2022-SupTab7-inversions-v1.bed", header = F)
colnames(Venu_2022_Inv) <- c("chr", "start", "end", "name")
Venu_2022_Inv$chr <- gsub("chr", "", Venu_2022_Inv$chr)

Venu_2022_Inv_v5 <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Venu-et-al-2022-SupTab7-inversions-v5.bed", header = F)
colnames(Venu_2022_Inv_v5) <- c("chr", "start", "end", "name","segment")
Venu_2022_Inv_v5$chr <- gsub("chr", "", Venu_2022_Inv_v5$chr)
Venu_2022_Inv_v5$Ecotype <- "NA"
Venu_2022_Inv_v5$Population <- "NA"


# Filer dataset to specific region
pca.comp.df.filt <- pca.comp.df %>%
  rowwise() %>%
  filter(any(
    chr == Venu_2022_Inv_v5$chr &
      start >= Venu_2022_Inv_v5$start-1e6 &
      end <= Venu_2022_Inv_v5$end+1e6
  )) %>%
  ungroup()

## Line graph of regions from Venu
regions_plot <- ggplot(pca.comp.df.filt,
                       aes(as.numeric(end), MDS1_scaled, col = Population)) +
  geom_segment(data = Venu_2022_Inv_v5, aes(x = start, xend = end, y = max(pca.comp.df.filt$MDS1_scaled), yend = max(pca.comp.df.filt$MDS1_scaled)),
                                         col = "red") +
  geom_line(aes(group = sample)) +
  facet_grid(.~chr,scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows_Venu_et_al_2022.png"), regions_plot, width = 40, height = 15)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows_Venu_et_al_2022.pdf"), regions_plot, width = 40, height = 15)

## Plot tile of Venu et al 2022 inversions
regions_plot <- ggplot(pca.comp.df.filt) +
  geom_segment(data = Venu_2022_Inv_v5, aes(x = as.numeric(start), xend = as.numeric(end), y = "Venu_et_al_2022", yend = "Venu_et_al_2022"),
                                         col = "red", size = 5) +       
  geom_tile(aes(as.numeric(end), sample, fill = MDS1_ratio)) +
  scale_fill_gradient2(low = "deepskyblue", mid = "orange" ,high = "darkgreen", midpoint=0.5) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5),
        panel.background = element_rect(fill = NA, color = "black"), 
        text = element_text(size = 20))

ggsave(paste0(plot.dir, pca_mds_file, "_mds_tile_specificWindows_Venu_et_al_2022.png"), regions_plot, width = 40, height = 15)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_tile_specificWindows_Venu_et_al_2022.pdf"), regions_plot, width = 40, height = 15)