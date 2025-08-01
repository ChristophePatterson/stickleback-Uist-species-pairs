source("~/apps/twisst/plot_twisst.R")
# Read in unique combination of populations
args <- commandArgs(trailingOnly = TRUE)
# set path
top_dir <- args[1]
# top_dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/ploidy_aware_HWEPops_MQ10_BQ20/twisst/Population_comparison"

setwd(top_dir)
library(ggplot2)
# devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)
library(patchwork)
library(tidyverse)

pop_combn <- read.table("pop_uniq.txt_combn.txt")
pop_combn$comb <- paste0(pop_combn[,1], "_", pop_combn[,2])
# pop_combn <- pop_combn[-grep("OLAV", pop_combn$comb),]

twisst_data_all <- list()
for(i in 1:length(pop_combn$comb)){
  run_name <- pop_combn$comb[i]

  #weights file with a column for each topology
  weights_file <- paste0(run_name,"/", run_name, ".weights.tsv.gz")
  
  # It is not necessary to import window data files, but if you do there should be one for
  # each weights file
  
  #coordinates file for each window
  window_data_file <- paste0(run_name,"/", run_name, ".data.tsv")
  
  ################################# import data ##################################
  
  # The function import.twisst reads the weights, window data  files into a list object
  # If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
  # in the window data file, these will be separated when importing.
  
  twisst_data <- import.twisst(weights_files=weights_file,
                                    window_data_files=window_data_file)

  #subset twisst object for these
  window_data_df <- do.call("rbind", twisst_data[["window_data"]])
  weight_data_df <- do.call("rbind", twisst_data[["weights_raw"]])
  window_data_df$run_name <- run_name
  weight_data_df$run_name
  twisst_data_all[[i]] <- cbind.data.frame(window_data_df, weight_data_df)
}

# Combine into single dataset
twisst_data_all <- do.call("rbind.data.frame", twisst_data_all)
# Colour blind pallette
cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")

## Swap out seq number for scaffold
scaf <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
twisst_data_all$chr <- scaf$Sequence.name[match(twisst_data_all$scaffold,scaf$RefSeq.seq.accession)]

write.table(twisst_data_all, file = paste0(top_dir, "/twisst_all_population_combinations.txt"))

p <- ggplot(twisst_data_all) +
  geom_step(aes(start, topo2, col = run_name), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x") +
  # scale_color_manual(values = topo_cols) +
  theme_classic()

ggsave(file = paste0(top_dir, "/All_population_comparison_twisst.pdf"), p, width = 15, height = 30)
ggsave(file = paste0(top_dir, "/All_population_comparison_twisst.png"), p, width = 15, height = 30)

## Create tree files
tree.text <- c("((Resident_A, Resident_B),(Migratory_A, Migratory_B));",
               "((Resident_A, Migratory_A),(Resident_B, Migratory_B));",
               "((Migratory_A, Resident_B),(Resident_A, Migratory_B));")

## Read trees
trees <- read.tree(text = tree.text)

# Set ecotype colour
eco.cols <- c("orange", "deepskyblue")

## Make ecotype tree red
p1 <- ggtree(trees[[1]], layout = "slanted", size = 2, col = "firebrick1") +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(0, 6) +
  scale_color_manual(values = eco.cols) +
  ggtitle("(a) Ecotype")
## Geographic tree
p2 <- ggtree(trees[[2]], layout = "slanted", size = 2) +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(0, 6)  +
  scale_color_manual(values = eco.cols) +
    ggtitle("(b) Geographic")

# Alternate tree
p3 <- ggtree(trees[[3]], layout = "slanted", size = 2) +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(0, 6)  +
  scale_color_manual(values = eco.cols) +
  ggtitle("(c) Alternate")

tree.plot <- p1 + p2 + p3 

ggsave(filename = "tree_comb.png", tree.plot, width = 11.5, height = 5)

# Remove chr from chr names
twisst_data_all$chr <- gsub("chr", "", twisst_data_all$chr)

## Create heat map for ecotype tree
p <- ggplot(twisst_data_all) +
  geom_tile(aes(mid, run_name, col = topo2)) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "white", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Ecotype Tree Weight") +
  facet_grid(chr~.) +# , scale = "free_x", space = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom") +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) 
  
# Combine with tree plot
twisst_tree_plot <- (tree.plot / p) + plot_layout(heights = c(1,10))

ggsave(filename = "twisst_combined.png", twisst_tree_plot, width = 10, height = 20)