# Set library paths
.libPaths(c("/gpfs01/home/mbzcp2/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

# Load libraries
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)

## Load in arguments
args <- commandArgs(trailingOnly = TRUE)
RAxML.run <- args[1]
# Get support tree name
tree.name <- paste0(RAxML.run, ".support")

## Get colour palette
# Colour blind friendly palette
cbPalette <- c("#D55E00","#009E73","#E69F00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# Read in tree
tree <- read.tree(tree.name)

# Read in sample metadata
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = TRUE)

## Plot tree
plot.tree <- ggtree(tree, aes(col = Waterbody), "ape") #"daylight")

# Merge metadata with tree
plot.tree <- plot.tree %<+% samples

# Recode ecotype for plotting
plot.tree$data$Ecotype[plot.tree$data$Ecotype=="anad"] <- "mig"
plot.tree$data$Ecotype[!plot.tree$data$Ecotype%in%c("mig","resi","st", "fw")] <- "xother"

# Create support value
plot.tree$data$support <- as.numeric(plot.tree$data$label)
# Create binary 100 percent bootstrap tree
plot.tree$data$support.100 <- plot.tree$data$support
plot.tree$data$support.100[plot.tree$data$support<95] <- NA
plot.tree$data$support.100[plot.tree$data$support>=95] <- "*"

## Custom tip colours
plot.tree <- plot.tree + geom_tippoint(aes(shape = Ecotype, col = Waterbody), size = 5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = cbPalette) +
  geom_nodepoint(aes(fill = support), size = 3, shape = 21, col = "black") +
  scale_fill_gradient(low = "grey",high = "firebrick3", name = "Bootstrap\nSupport (%)") +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "firebrick3", size = 10) 

# ggsave("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/test.png", plot.tree, width = 7.96, height = 7.96)
ggsave(paste0(RAxML.run, ".plot.png"), plot.tree, width = 7.96, height = 7.96)
ggsave(paste0(RAxML.run, ".plot.pdf"), plot.tree, width = 7.96, height = 7.96)
