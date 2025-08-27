source("~/apps/twisst/plot_twisst.R")
# Read in unique combination of populations
args <- commandArgs(trailingOnly = TRUE)
# set path
# top_dir <- args[1]
top_dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/twisst/Population_comparison/"

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
  weights_file <- paste0(run_name,"/", run_name, ".all.weights.tsv.gz")
  
  # It is not necessary to import window data files, but if you do there should be one for
  # each weights file
  
  #coordinates file for each window
  window_data_file <- paste0(run_name,"/", run_name, ".all.data.tsv")
  
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

## Remove data that equals NA
table(apply(twisst_data_all, MARGIN = 1, function(x) max.topo <- any(is.na(x))))
twisst_data_all <- twisst_data_all[!apply(twisst_data_all, MARGIN = 1, function(x) max.topo <- any(is.na(x))), ]
# Colour blind pallette
cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")

## Swap out seq number for scaffold
scaf <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
# Remove small contigs and mitogenome
scaf <- scaf[!scaf$Chromosome.name%in%c("Un","MT"),]
# Get chromosome name
twisst_data_all$chr <- scaf$Sequence.name[match(twisst_data_all$scaffold, scaf$GenBank.seq.accession)]

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

## Calculate proportion of weighting across all species pairs

twisst_data_all$top_weight_topo <- apply(twisst_data_all, MARGIN = 1, function(x) max.topo <- c("topo1","topo2", "topo3")[which.max(x[c("topo1","topo2", "topo3")])])

topo_weight_counts <- twisst_data_all %>%
  group_by(run_name) %>%
  group_by(top_weight_topo, .add = T) %>%
  summarise(count = n()) %>%
  group_by(run_name) %>%
  mutate(freq = count/sum(count)) %>%
  group_by(top_weight_topo) %>%
  summarise(mn.topo = mean(freq, na.rm=T))

tree_weight_sum <- twisst_data_all %>%
  group_by(run_name) %>%
  summarise(mn_topo1 = mean(topo1, na.rm = T),
            mn_topo2 = mean(topo2, na.rm = T),
            mn_topo3 = mean(topo3, na.rm = T),
            md_topo1 = median(topo1, na.rm = T),
            md_topo2 = median(topo2, na.rm = T),
            md_topo3 = median(topo3, na.rm = T),
            mn_dist_topo1v2 = mean(topo2-topo1, na.rm = T),
            mn_dist_topo1v3 = mean(topo2-topo3, na.rm = T),
            top_topo1.perc = (sum(top_weight_topo=="topo1")/length(top_weight_topo))*100,
            top_topo2.perc = (sum(top_weight_topo=="topo2")/length(top_weight_topo))*100,
            top_topo3.perc = (sum(top_weight_topo=="topo3")/length(top_weight_topo))*100)

tree_weight_greater_perc <- twisst_data_all %>%
  group_by(run_name) %>%
  summarise(topo1_gr_0.5 = sum(topo1>=(1/2))/length(topo1)*100,
            topo2_gr_0.5 = sum(topo2>=(1/2))/length(topo2)*100,
            topo3_gr_0.5 = sum(topo3>=(1/2))/length(topo3)*100,
            topo1_gr_0.666 = sum(topo1>=(2/3))/length(topo1)*100,
            topo2_gr_0.666 = sum(topo2>=(2/3))/length(topo2)*100,
            topo3_gr_0.666 = sum(topo3>=(2/3))/length(topo3)*100)

tree_weight_greater_perc
summary(tree_weight_greater_perc)

## Which regions of the genome are favoured greater than 2/3 for all population comparison
twisst_data_all$wndname <- paste0(twisst_data_all$scaffold, "_", twisst_data_all$start)

top_weight_all_comp <- twisst_data_all %>%
    group_by(wndname) %>%
    summarise(
              start = first(start),
              end = first(end),
              mid = first(mid),
              sites = first(sites),
              chr = first(chr),
              scaffold = first(scaffold),
              min.topo1 = min(topo1, na.rm = T),
              min.topo2 = min(topo2, na.rm = T),
              min.topo3 = min(topo3, na.rm = T),
              all_eco_over_half = all(topo2>=1/2),
              all_eco_over_2.3rds = all(topo2>=2/3),
              all_eco_over_95 = all(topo2>=95/100),
              all_geo_over_half = all(topo1>=1/2),
              all_alt_over_half = all(topo3>=1/2),
              all_geo_over_2.3rds = all(topo1>=2/3),
              all_alt_over_2.3rds = all(topo3>=2/3),
              all_geo_over_95 = all(topo1>=95/100),
              all_alt_over_95 = all(topo3>=95/100)
              )

top_weight_all_comp %>%
  summarise(all_eco_over_half = sum(all_eco_over_half)/nrow(top_weight_all_comp)*100,
            all_eco_over_2.3rds =sum(all_eco_over_2.3rds)/nrow(top_weight_all_comp)*100,
            all_eco_over_95 = sum(all_eco_over_95)/nrow(top_weight_all_comp)*100,
            all_geo_over_half = sum(all_geo_over_half)/nrow(top_weight_all_comp)*100,
            all_alt_over_half = sum(all_alt_over_half)/nrow(top_weight_all_comp)*100,
            all_geo_over_2.3rds = sum(all_geo_over_2.3rds)/nrow(top_weight_all_comp)*100,
            all_alt_over_2.3rds = sum(all_alt_over_2.3rds)/nrow(top_weight_all_comp)*100,
            all_geo_over_95 = sum(all_geo_over_95)/nrow(top_weight_all_comp)*100,
            all_alt_over_95 = sum(all_alt_over_95)/nrow(top_weight_all_comp)*100,
  ) %>%
  t()

max(top_weight_all_comp$min.topo1)
max(top_weight_all_comp$min.topo3)

twisst_data_all %>% 
  group_by(run_name) %>%
  summarise(mx.topo1 = max(topo1, na.rm = T),
            mx.topo2 = max(topo2, na.rm = T),
            mx.topo3 = max(topo3, na.rm = T))

# Create barplot
pbar <- ggplot(twisst_data_all) +
  geom_bar(aes(gsub("_", "-", run_name), fill = factor(top_weight_topo, levels = c("topo3","topo1", "topo2"))), position = "fill", show.legend = F) +
  scale_fill_manual(values = c("grey60","firebrick1","grey40"), breaks = c("topo1","topo2", "topo3"), labels = c("Geo","Eco", "Alt")) +
  coord_flip() +
  labs(fill = "Tree", x = "Waterbody Pairs", y = "") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) + 
  theme(panel.grid = element_blank(), legend.position = "bottom", axis.text.y=element_text(angle=-45, vjust = 1, hjust=1))

ggsave(filename = "Top_tree_topo_all_pop_combs.png", pbar, width = 5, height = 5)

## Make ecotype tree red
p1 <- ggtree(trees[[1]], layout = "slanted", size = 2, col = "firebrick1") +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(-2, 6) +
  scale_color_manual(values = eco.cols) +
  ggtitle("Ecotype")
## Geographic tree
p2 <- ggtree(trees[[2]], layout = "slanted", size = 2, col = "grey60") +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(-2, 6)  +
  scale_color_manual(values = eco.cols) +
    ggtitle("Geographic")

# Alternate tree
p3 <- ggtree(trees[[3]], layout = "slanted", size = 2, col = "grey40") +
  geom_tiplab(aes(col = substr(label, 1, 3)), show.legend = F, hjust = -0.1) +
  xlim(-2, 6)  +
  scale_color_manual(values = eco.cols) +
  ggtitle("Alternate")

tree.plot <- p1 + p2 + p3 + plot_annotation()
ggsave(filename = "tree_comb.png", tree.plot, width = 11.5, height = 5)

p.hist <- ggplot(twisst_data_all) +
  geom_histogram(aes(topo1), , alpha = 0.5, fill = "grey60", bins = 100) +
  geom_histogram(aes(topo2), , alpha = 0.5, fill = "firebrick1", bins = 100) +
  geom_histogram(aes(topo3), , alpha = 0.5, fill = "grey40", bins = 100) +
  theme_bw() 

p.hist.z <- ggplot(twisst_data_all) +
  geom_histogram(aes(topo1), , alpha = 0.5, fill = "grey60", bins = 100) +
  geom_histogram(aes(topo2), , alpha = 0.5, fill = "firebrick1", bins = 100) +
  geom_histogram(aes(topo3), , alpha = 0.5, fill = "grey40", bins = 100) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 100), clip="on") +
  facet_wrap(~run_name)


ggsave(filename = "Top_tree_topo_histo.png", p.hist/p.hist.z, width = 8, height = 8)

# Remove chr from chr names
twisst_data_all$chr <- gsub("chr", "", twisst_data_all$chr)
top_weight_all_comp$chr <- gsub("chr", "", top_weight_all_comp$chr)
twisst_data_all$chr <- factor(twisst_data_all$chr, levels = gsub("chr", "", scaf$Sequence.name[order(scaf$GenBank.seq.accession)]))
top_weight_all_comp$chr <- factor(top_weight_all_comp$chr, levels = gsub("chr", "", scaf$Sequence.name[order(scaf$GenBank.seq.accession)]))

## Create heat map for ecotype tree
pEco <- ggplot(twisst_data_all) +
  geom_segment(data = top_weight_all_comp[top_weight_all_comp$all_eco_over_2.3rds,], aes(x = start, xend = end, " "), col = "firebrick3", linewidth = 1, lineend = "round") +
  geom_segment(data = top_weight_all_comp[top_weight_all_comp$all_eco_over_95,], aes(x = start, xend = end, " "), col = "black", linewidth = 1, lineend = "round") +
  geom_segment(aes(x = start, xend = end, run_name, col = topo2), linewidth = 2.8) +
  coord_cartesian(ylim = c(1, 6), clip="off") +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick3", midpoint = 1/3, name = "Ecotype Tree Weight", limits = c(0, 1)) +
  facet_grid(chr~.) +# , scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey", color = "black", linewidth = 0.25), legend.position = "bottom",
        axis.text.y = element_text(size = 6), axis.ticks.y = element_blank(), strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
        axis.line = element_blank()) +
  ylab("Waterbody Pairs") 


# Combine with tree plot
twisst_tree_plot <- (tree.plot / pbar / pEco ) + plot_layout(heights = c(2, 1, 10), tag_level = "new") + plot_annotation(tag_levels = list(c('(a)', '', '', "(b)","(c)")))
ggsave(filename = "twisst_combined.png", twisst_tree_plot , width = 10, height = 20)
# ggsave(filename = "/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/test.png", twisst_tree_plot , width = 10, height = 20)


## Create heat map for ecotype tree
pEco <- ggplot(twisst_data_all) +
  geom_segment(data = top_weight_all_comp, aes(x = start, xend = end, chr, col = min.topo2), linewidth = 7) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Ecotype Tree Weight", limits = c(0, 1)) +
  # facet_grid(chr~.) +# , scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom",
        axis.text.y = element_text(size = 6)) +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  ylab("Waterbody Pairs") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) 
  
# Combine with tree plot
twisst_tree_plot <- (tree.plot / pbar / pEco ) + plot_layout(heights = c(2, 1, 5), tag_level = "new") + plot_annotation(tag_levels = list(c('(a)', '', '', "(b)","(c)")))

ggsave(filename = "twisst_combined_min_single track.png", twisst_tree_plot , width = 10, height = 20)


# Geographic heat map
pGeo <- ggplot(twisst_data_all) +
  geom_segment(aes(x = start, xend = end, run_name, col = topo1), linewidth = 4) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Geographic Tree Weight", limits = c(0, 1)) +
  facet_grid(chr~.) +# , scale = "free_x", space = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom", axis.text.y = element_blank(),axis.title.y = element_blank()) +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  ylab("Waterbody Pair") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) 

## Alternate heatmap
pAlt <- ggplot(twisst_data_all) +
  geom_segment(aes(x = start, xend = end, run_name, col = topo3), linewidth = 4) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Alternate Tree Weight", limits = c(0, 1)) +
  facet_grid(chr~.) +# , scale = "free_x", space = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom", axis.text.y = element_blank(),axis.title.y = element_blank()) +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  ylab("Waterbody Pair") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) 

## Combine with tree plots
pall <- ((p1 + p2 + p3) / (pEco + pGeo + pAlt)) + plot_layout(heights =c(1,10))

# Save
ggsave(filename = "twisst_combined_alltopos.png", pall, width = 20, height = 20)

## Zoom in on regions of interest
regions <- data.frame(chr = c("I", "IX", "XI", "XXI"), start = c(25000000, 4500000, 5000000, 8000000), end = c(31000000, 10000000, 10000000, 15000000))

# Filer dataset to specific region
twisst_data_all_filt <- twisst_data_all %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

## Create heat map for ecotype tree for regions of interest
pEco_zoom <- ggplot(twisst_data_all_filt) +
  geom_segment(aes(x = start, xend = end, run_name, col = topo2), linewidth = 12) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Ecotype Tree Weight", limits = c(0, 1)) +
  facet_wrap(~chr, ncol = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom") +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  ylab("Waterbody Pair") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(twisst_data_all$end),1e6)),name = "Mbs",expand = expansion(0)) 
  
# Combine with tree plot
twisst_tree_plot <- (tree.plot / pEco_zoom) + plot_layout(heights = c(1,6))

ggsave(filename = "twisst_combined_zoom.png", twisst_tree_plot, width = 10, height = 20)

## Select just chri
chrI_select <- data.frame(chr = c("I"), start = c(26500000), end = c(27200000))

# Filer dataset to specific region
twisst_data_all_filt_chrI <- twisst_data_all %>%
  rowwise() %>%
  filter(any(
    chr == chrI_select$chr &
      start >= chrI_select$start &
      end <= chrI_select$end
  )) %>%
  ungroup()

twisst_data_all_filt_chrI

pEco_zoom_chrI <- ggplot(twisst_data_all_filt_chrI) +
  geom_segment(aes(x = start, xend = end, run_name, col = topo2), linewidth = 12) +
  # scale_color_viridis_c(option = "rocket") +
  scale_color_gradient2(low = "black", mid =  "white", high = "firebrick1", midpoint = 1/3, name = "Ecotype Tree Weight", limits = c(0, 1)) +
  # facet_wrap(~chr, ncol = 1, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey"), legend.position = "bottom") +
  guides(colour = guide_colorbar(theme = theme(legend.frame = element_rect(colour = "black")))) +
  ylab("Waterbody Pairs") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(chrI_select$start, chrI_select$end, 1e5)),name = "Mbs",expand = expansion(0)) 


ggsave(filename = "twisst_combined_zoom_chrI.png", pEco_zoom_chrI, width = 6, height = 4)
