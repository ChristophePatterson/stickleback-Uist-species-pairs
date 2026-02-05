library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly=T)

plot.dir <- args[1]
# plot.dir <- "C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/pca/Anad_resi/wndsize25000_wndslid5000/"
# plot.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/pca/Anad_resi_fw/wndsize25000_wndslid5000/"

pca_mds_file <- gsub(".txt", "", args[2])
# pca_mds_file <- gsub(".txt", "", "sliding-window_pca_wndsize25000_wndslid5000.txt")

pca.comp.df <- read_csv(paste0(plot.dir, pca_mds_file, ".txt"))
# pca.comp.df <- read_csv("C:/Users/mbzcp2/Downloads/Anad_resi/wndsize25000_wndslid5000/sliding-window_pca_wndsize25000_wndslid5000.txt")

# Read in chromosome details
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
# Remove smaller scafs and mito
chr <- chr[!chr$Chromosome.name%in%c("Un","MT"),]

chr$Sequence.name <- gsub("chr", "", chr$Chromosome.name) %>%
  as.roman() %>%
  as.character()
# Order Chromsome name to be in order 
chr$Sequence.name <- factor(chr$Sequence.name, levels = chr$Sequence.name[order(chr$GenBank.seq.accession)])

# Calculate cumulative position of start of chromosome
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))

chr$GenBank.seq.accession[1:10]

# Order chr in order in pca
pca.comp.df$chr <- chr$Sequence.name[match(pca.comp.df$chr, chr$GenBank.seq.accession)]

pca.comp.df$chr[1:10]
pca.comp.df$chr <- factor(pca.comp.df$chr, levels = chr$Sequence.name[order(chr$GenBank.seq.accession)])

# Create column each window for eaiser subset
pca.comp.df$windowname <- paste(pca.comp.df$chr, pca.comp.df$start, pca.comp.df$end, sep = "-")

## Add in sample data
sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
pca.comp.df$Population <- sample_data$Population[match(pca.comp.df$sample, sample_data$individual)]
pca.comp.df$Ecotype <- sample_data$Ecotype[match(pca.comp.df$sample, sample_data$individual)]
pca.comp.df$Waterbody <- sample_data$Waterbody[match(pca.comp.df$sample, sample_data$individual)]

# Remove stream
pca.comp.df$Ecotype[pca.comp.df$Ecotype=="st"] <- "fw"
pca.comp.df$Ecotype[pca.comp.df$Ecotype=="anad"] <- "mig"
## Change level
pca.comp.df$Ecotype <- factor(pca.comp.df$Ecotype, levels = c("mig", "resi", "fw")) 

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

# Compute the sign for each window and ecotype == "mig"
signs <- pca.comp.df %>%
  filter(Ecotype == "mig") %>%
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

#### ## Plot PCA 1 across the genome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
####   geom_point() +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## Plot PCA 2 across the genome
#### q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
####   geom_point() +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### # Save output
#### ## ggsave(paste0(plot.dir, pca_mds_file,"_pca12_point.png"), p/q, width = 40, height = 15)
#### 
#### ## Plot PCA 1 across the genome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## Plot PCA 2 across the genome
#### q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### # Save output
#### ## ggsave(paste0(plot.dir, pca_mds_file,"_pca12_line.png"), p/q, width = 40, height = 15)
#### 
#### # Break down by chromosome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
####   geom_point() +
####   facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)), name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file,"_Popsplit.png"), p, width = 40, height = 20)
#### 
#### ## PLot MDS1 along genome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
####   geom_point() +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds1.png"), p, width = 40, height = 20)
#### 
#### ## PLot MDS1 along genome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(chr~., scale = "free_x", space = "free_x", switch = "y") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds1_vert.png"), p, width = 30, height = 40)
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds1_vert.pdf"), p, width = 30, height = 40)
#### 
#### ## Plot MDS1 along genome and split by waterbody
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
####   geom_point() +
####   facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "mds_Popsplit.png"),p , width = 40, height = 20)
#### 
#### ## Plot MDS1 along genome and split by waterbody
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "line_Popsplit.png"), p, width = 40, height = 20)
#### ## ggsave(paste0(plot.dir, pca_mds_file, "line_Popsplit.pdf"), p, width = 40, height = 20)
#### 
#### ## Plot PCA 1 across the genome
#### p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## Plot PCA 2 across the genome
#### q <- ggplot(pca.comp.df, aes(as.numeric(end), MDS2_scaled, col = Population, shape = Ecotype)) +
####   geom_line(aes(group = sample)) +
####   facet_grid(.~chr, scale = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### # Save output
#### ## ggsave(paste0(plot.dir, pca_mds_file,"mds_line.png"), p/q, width = 40, height = 15)


### Zoomed in sections
# regions <- data.frame(chr = c("I", "IX", "XI", "XXI"), start = c(25000000, 4500000, 5000000, 8000000), end = c(31000000, 10000000, 10000000, 15000000))
regions <- data.frame(chr = c("I", "IX", "XI", "XXI"), start = c(26300000, 5500000, 6000000, 9400000), end = c(27500000, 8700000, 7000000, 11500000))
regions$start.cum <- regions$start+(chr$Cum.Seq.length[match(regions$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])


# Filer dataset to specific region
pca.comp.df.filt <- pca.comp.df %>%
  filter(chr %in% regions$chr) %>%
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

## ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows.png"), regions_plot, width = 40, height = 15)
## ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows.pdf"), regions_plot, width = 40, height = 15)

## Segement plot  of regions
genome.all.plot <- ggplot(chr) +
  geom_segment(aes(x = Cum.Seq.length, xend = Cum.Seq.length+Seq.length, y = 0, col = as.numeric(Chromosome.name) %% 2 == 0), linewidth = 500) +
  geom_segment(data = regions, aes(x = start.cum, xend = start.cum+(start-end), y = 0), linewidth = 500, col = "firebrick") +
  geom_text(aes(Cum.Seq.length+(Seq.length/2), y = 0, label = Sequence.name), col = "white") + 
  geom_vline(aes(xintercept = Cum.Seq.length)) +
  scale_color_manual(values = c("grey30", "grey50")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), label = function(x) x/1e6, name = "Mbps") +
  theme_void() +
  theme(legend.position = "none")
  #theme(legend.position = "none", axis.text.x = element_text(), axis.title.x = element_text(),
  #      plot.margin = rep(unit(10, units = "points"), 4))
genome.all.plot


tile_plot <- ggplot(pca.comp.df.filt[pca.comp.df.filt$chr!="I",],
                    aes(as.numeric(end), sample, fill = MDS1_ratio, shape = Ecotype)) +
  geom_tile() +
  #scale_fill_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", midpoint=0.5, name =  "MDS Scaled") +
  scale_fill_gradient2(low = "#009E73", mid = "#E69F00", high = "#56B4E9", midpoint=0.5, name =  "MDS Scaled") +
  #scale_fill_gradient2(low = "firebrick3", mid = "orange" ,high = "darkgreen", midpoint=0.5, name =  "MDS Scaled") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),0.5e6)),name = "Mbps", expand = c(0.01,0)) +
  facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "top", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(0.5,'lines'),
        legend.frame = element_rect(colour="black"),
        legend.ticks = element_line(colour="black"),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
        axis.text.y = element_blank(),                    # hide left/right axis ticks
        axis.title.y = element_blank(), 
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(fill = "grey90", color = "black", size = 0.5))

tile_plot_full <- tile_plot / genome.all.plot + plot_layout(heights = c(10,1))

ggsave(paste0(pca_mds_file, "_mds_ratio_tile_specificWindows.png"), tile_plot_full, height = 7.96*1.1, width = 24.62*0.66666)

ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_specificWindows.png"), tile_plot_full, height = 7.96*1.1, width = 24.62*0.66666)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_specificWindows.pdf"), tile_plot_full, height = 7.96*1.1, width = 24.62*0.66666)

## Location of ATP1A1
ATP1A1 <- data.frame(chr = "I", start = 26848100, end = 26861753, name = "atp1a1a", Ecotype="", Population="")

# Plot zoomed inversion on chrI
tile_plot_chrI <- ggplot(pca.comp.df[pca.comp.df$chr=="I"&pca.comp.df$start>=26500000&pca.comp.df$end<=27200000,]) +
  geom_tile(aes(as.numeric(end), sample, fill = MDS1_ratio)) +
  geom_vline(xintercept = as.numeric(ATP1A1$start), col = "black", size = 0.5) + 
  geom_vline(xintercept = as.numeric(ATP1A1$end), col = "black", size = 0.5) +  
  # add gene label on bottom of plot outside of tiles
  annotate("text", x = (as.numeric(ATP1A1$start) + as.numeric(ATP1A1$end))/2, y = 1, label = ATP1A1$name, color = "black", size = 4, vjust = 1.5) +
  # scale_fill_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", midpoint=0.5, name =  "MDS Scaled") +
  scale_fill_gradient2(low = "#009E73", mid = "#E69F00", high = "#56B4E9", midpoint=0.5, name =  "MDS Scaled") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),0.1e6)),name = "Mbps", expand = c(0.01,0)) +
  facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  coord_cartesian(clip = "off") +  # allow drawing in the margin area
  theme(legend.position = "top", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(10,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
        axis.text.y = element_blank(),                    # hide left/right axis ticks
        axis.title.y = element_blank(), 
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(fill = "grey90", color = "black", size = 0.5))


tile_plot_full <- genome.all.plot / tile_plot_chrI / tile_plot + guide_area() + plot_layout(heights = c(1, 10, 10,1), guides = "collect")

ggsave(paste0("test.png"), tile_plot_full , height = 7.96*1.5, width = 24.62/2)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_chI_inv.png"), tile_plot_chrI , height = 7.96*1.1, width = 24.62*0.66666)
ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_chI_inv.pdf"), tile_plot_chrI , height = 7.96*1.1, width = 24.62*0.66666)


# Filer dataset to specific region
pca.comp.df.filt.within <- pca.comp.df %>%
  filter(chr %in% regions$chr) %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start+(regions$end-regions$start)*0.3 &
      end <= regions$end-(regions$end-regions$start)*0.3
  )) %>%
  ungroup()


haplotype_thresholds <- c(0.33, 0.66)
hapltotypes <- pca.comp.df.filt.within %>%
  group_by(chr, sample) %>%
  summarise(median_mds1 = median(MDS1_ratio)) %>%
  mutate(haplotype = case_when(
    median_mds1 < haplotype_thresholds[1] ~ "Derived",
    median_mds1 >= haplotype_thresholds[1] & median_mds1 < haplotype_thresholds[2] ~ "Heterozygous",
    median_mds1 >= haplotype_thresholds[2] ~ "Ancestral"),
    individual = sample) %>%
  left_join(sample_data) %>%
  mutate(Ecotype = case_when(Ecotype == "anad" ~ "Mig",
                             Ecotype == "resi" ~ "Resi"),
         haplotype = factor(haplotype, levels = c("Derived", "Heterozygous", "Ancestral")))

hist(hapltotypes$median_mds1)

hap.plot <- ggplot(hapltotypes) +
  geom_bar(aes(y = Ecotype, fill = haplotype), position = "stack") +
  facet_grid(Waterbody~chr, scales = "free_y") +
  theme_classic() +
  #scale_x_continuous(expand = c(0.1,0), position = "right") +
  scale_fill_manual(values = c("#009E73", "#E69F00","#56B4E9"), name = "Haplotype") + 
  theme(legend.position = "bottom", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(1,'lines'),
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(color = "black", size = 0.5))
  
tile_plot_full <- genome.all.plot + 
  (tile_plot_chrI  + theme(legend.position = "none")) + 
  (tile_plot + theme(legend.position = "bottom")) + 
  (hap.plot + theme(legend.position = "bottom")) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(heights = c(1,10,10),widths = c(3,2), 
  design = c("
              AD
              BD
              CC
              "))

ggsave(paste0(plot.dir, pca_mds_file, "_full_inv_plot.png"), tile_plot_full , height = 7.96*1.6, width = 24.62*0.75)
ggsave(paste0(plot.dir, pca_mds_file, "_full_inv_plot.png"), tile_plot_full , height = 7.96*1.6, width = 24.62*0.75)

### Viz of pca1 and pca2 surronding atp1a1
atp1a1_region <- pca.comp.df %>%
  filter(chr == ATP1A1$chr & start >= (ATP1A1$start-50000) & end <= (ATP1A1$end+50000))

## number of windows
nrow(atp1a1_region[atp1a1_region$sample==atp1a1_region$sample[1],])

## plot each window separately
regions_plots <- list()
for(i in 1:nrow(regions)){
  window_data <- regions[i, ]
  window_pca <- pca.comp.df %>%
    filter(chr == window_data$chr & start >= window_data$start & end <= window_data$end)
  regions_plots[[i]] <- ggplot(window_pca, aes(start, sample, fill = MDS1_ratio)) +
    geom_tile(size = 3) +
    theme_bw() +
    scale_fill_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", midpoint=0.5, name =  "MDS Scaled") +
    ggtitle(paste0("Window: ", window_data$chr[1], ":", window_data$start[1], "-", window_data$end[1])) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),0.1e6)),name = "Mbs", expand = c(0.01,0)) +
    facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
    theme_classic() +
    coord_cartesian(clip = "off") +  # allow drawing in the margin area
    theme(legend.position = "bottom", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(0.5,'lines'),
          axis.title.y.right = element_blank(),                # hide right axis title
          axis.text.y.right = element_blank(),                 # hide right axis labels
          axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
          axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
          strip.background = element_rect(color = "black", size = 0.5),
          panel.background = element_rect(fill = "grey90", color = "black", size = 0.5))
  ## Save plot
  ggsave(paste0(plot.dir, pca_mds_file, "_mds_ratio_tile_", window_data$chr[1], "_wnd", window_data$start[1]/1000, "kb_sld", window_data$end[1]/1000, "kb.png"), regions_plots[[i]], width = 15, height = 10)
  
}

## Calculate which haplotype each sample belongs to for chrIX inversion
# Get median MDS1 value for each sample in the inversion region
inversion_chrIX_region <- pca.comp.df %>%
  filter(chr == "IX" & start >= 5900000 & end <= 7400000)
median_mds1_by_sample <- inversion_chrIX_region %>%
  group_by(sample) %>%
  summarise(median_mds1 = median(MDS1_ratio, na.rm = TRUE)) 

# Define haplotype thresholds
haplotype_thresholds <- c(0.33, 0.66)
# Classify samples into haplotypes
median_mds1_by_sample <- median_mds1_by_sample %>%
  mutate(haplotype = case_when(
    median_mds1 < haplotype_thresholds[1] ~ "AA",
    median_mds1 >= haplotype_thresholds[1] & median_mds1 < haplotype_thresholds[2] ~ "AT",
    median_mds1 >= haplotype_thresholds[2] ~ "TT"))

## Merge with sample data
median_mds1_by_sample <- median_mds1_by_sample %>%
  left_join(sample_data, by = c("sample" = "individual"))

## Merge with Genomic sex data
median_mds1_by_sample <- median_mds1_by_sample %>%
  left_join(read.table("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/Genomic_sex_determination.txt", header = F), by = c("sample" = "V1")) %>%
  rename(Genomic_sex = V2)

## Rename stream as freshwater
median_mds1_by_sample$Ecotype[median_mds1_by_sample$Ecotype=="st"] <- "fw"
median_mds1_by_sample$Ecotype[median_mds1_by_sample$Ecotype=="anad"] <- "mig"

# Plot histogram of median MDS1 values
hist_plot <- ggplot(median_mds1_by_sample, aes(x = median_mds1, fill = haplotype)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  facet_wrap(~Genomic_sex*Ecotype) +
  theme_bw() 

ggsave("test.png", hist_plot, width = 8, height = 6)

# Contingency tables
chrIX_Mig_haplotype_table <- table(median_mds1_by_sample$Genomic_sex[median_mds1_by_sample$Ecotype=="mig"], median_mds1_by_sample$haplotype[median_mds1_by_sample$Ecotype=="mig"])

# Chi-squared test between Males and females in migratory ecotype
chrIX_Mig_haplotype_chisq <- chisq.test(chrIX_Mig_haplotype_table)
# Print results
chrIX_Mig_haplotype_chisq  

#### ## Inversion Identified by Venu 2022
#### Venu_2022_Inv <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Venu-et-al-2022-SupTab7-inversions-v1.bed", header = F)
#### colnames(Venu_2022_Inv) <- c("chr", "start", "end", "name")
#### Venu_2022_Inv$chr <- gsub("chr", "", Venu_2022_Inv$chr)
#### 
#### Venu_2022_Inv_v5 <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Venu-et-al-2022-SupTab7-inversions-v5.bed", header = F)
#### colnames(Venu_2022_Inv_v5) <- c("chr", "start", "end", "name","segment")
#### Venu_2022_Inv_v5$chr <- gsub("chr", "", Venu_2022_Inv_v5$chr)
#### Venu_2022_Inv_v5$Ecotype <- "NA"
#### Venu_2022_Inv_v5$Population <- "NA"
#### 
#### 
#### # Filer dataset to specific region
#### pca.comp.df.filt <- pca.comp.df %>%
####   rowwise() %>%
####   filter(any(
####     chr == Venu_2022_Inv_v5$chr &
####       start >= Venu_2022_Inv_v5$start-1e6 &
####       end <= Venu_2022_Inv_v5$end+1e6
####   )) %>%
####   ungroup()
#### 
#### ## Line graph of regions from Venu
#### regions_plot <- ggplot(pca.comp.df.filt,
####                        aes(as.numeric(end), MDS1_scaled, col = Population)) +
####   geom_segment(data = Venu_2022_Inv_v5, aes(x = start, xend = end, y = max(pca.comp.df.filt$MDS1_scaled), yend = max(pca.comp.df.filt$MDS1_scaled)),
####                                          col = "red") +
####   geom_line(aes(group = sample)) +
####   facet_grid(.~chr,scales = "free_x", space = "free_x") +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows_Venu_et_al_2022.png"), regions_plot, width = 40, height = 15)
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds_line_specificWindows_Venu_et_al_2022.pdf"), regions_plot, width = 40, height = 15)
#### 
#### ## Plot tile of Venu et al 2022 inversions
#### regions_plot <- ggplot(pca.comp.df.filt) +
####   geom_segment(data = Venu_2022_Inv_v5, aes(x = as.numeric(start), xend = as.numeric(end), y = "Venu_et_al_2022", yend = "Venu_et_al_2022"),
####                                          col = "red", size = 5) +       
####   geom_tile(aes(as.numeric(end), sample, fill = MDS1_ratio)) +
####   scale_fill_gradient2(low = "deepskyblue", mid = "orange" ,high = "darkgreen", midpoint=0.5) +
####   scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
####   facet_grid(Ecotype+Population~chr,scale = "free", space = "free", switch = "y") +
####   theme_classic() +
####   theme(legend.position = "top",panel.spacing = unit(0,'lines'),
####         axis.title.y.right = element_blank(),                # hide right axis title
####         axis.text.y.right = element_blank(),                 # hide right axis labels
####         axis.ticks.y = element_blank(),                      # hide left/right axis ticks
####         axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
####         strip.background = element_rect(size = 0.5),
####         panel.background = element_rect(fill = NA, color = "black"), 
####         text = element_text(size = 20))
#### 
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds_tile_specificWindows_Venu_et_al_2022.png"), regions_plot, width = 40, height = 15)
#### ## ggsave(paste0(plot.dir, pca_mds_file, "_mds_tile_specificWindows_Venu_et_al_2022.pdf"), regions_plot, width = 40, height = 15)