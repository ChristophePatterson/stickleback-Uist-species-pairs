source("~/apps/twisst/plot_twisst.R")
library(tidyverse)

## Get commnad arguments
args <- commandArgs(trailingOnly = TRUE)

# Get data
run_name <- args[1]

#weights file with a column for each topology
weights_file <- paste0(run_name, ".weights.tsv.gz")

# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
window_data_file <- paste0(run_name, ".data.tsv")

################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


#################### subset to only the most abundant topologies #################

#get list of the most abundant topologies (top 4 in this case)
# Or if less than 4 use max number of topos

if(length(twisst_data$topos)<4){ topo.num <- length(twisst_data$topos)}
if(length(twisst_data$topos)>=4){ topo.num <- 4}
1:topo.num
top4_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:topo.num]

#subset twisst object for these
twisst_data_top4topos <- subset.twisst.by.topos(twisst_data, top4_topos)
#this can then be used in all the same plotting functions above.

# Convert twisst data into tidy format
window_data_df <- do.call("rbind", twisst_data_top4topos[["window_data"]])
weight_data_df <- do.call("rbind", twisst_data_top4topos[["weights_raw"]])
# Combine window and weight data
twisst_data_sf <- cbind(window_data_df, weight_data_df)

head(twisst_data_sf)
names(twisst_data_top4topos$topos)

# Which topo is favoured in each segment
twisst_data_sf$max.topo <- names(twisst_data_top4topos$topos)[max.col(twisst_data_sf[,names(twisst_data_top4topos$topos)])]
# Calclate runs of max topos
repeat.seq.topos <- rle(twisst_data_sf$max.topo)
twisst_data_sf$topo.max.repeats <- rep(repeat.seq.topos$lengths >= 5,times = repeat.seq.topos$lengths)
# Calculate max topo weighting
twisst_data_sf$max.topo.weight <- apply(twisst_data_sf[,names(twisst_data_top4topos$topos)], MARGIN = 1, function(x) max(x))
twisst_data_sf$topo_sig <- twisst_data_sf$topo.max.repeats&(twisst_data_sf$max.topo.weight>=0.5)
# Convert into tidy format
twisst_data_sf <- pivot_longer(twisst_data_sf, cols = names(twisst_data_top4topos$topos), names_to = "topo", values_to = "weight")
window_ave <- mean(twisst_data_sf$end - twisst_data_sf$start, na.rm = T)

### 
pdf(file = paste0(run_name, ".summary.pdf"), width  = 20, height = 40)
plot.twisst.summary(twisst_data_top4topos, lwd=3, cex=0.7)
dev.off()

pdf(file = paste0(run_name, ".all_typo.pdf"), width  = 20, height = 40)
plot.twisst(twisst_data_top4topos, mode=1, show_topos=TRUE)
dev.off()

pdf(file = paste0(run_name, ".all_typo_stacked.pdf"), width  = 40, height = 40)
plot.twisst(twisst_data_top4topos, mode=3, show_topos=TRUE)
dev.off()
### 
### # make smooth weightings and plot those across chromosomes
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
chr_select <- chr$GenBank.seq.accession[chr$Seq.length>10000000&!grepl("Y", chr$Chromosome.name)]
twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions = chr_select)

# Replace scaffold name with chr name
twisst_data_sf$chr <- chr$Chromosome.name[match(twisst_data_sf$scaffold, chr$GenBank.seq.accession)]

p <- ggplot(twisst_data_sf) +
  #geom_bar(stat="identity", position = "stack", width = 1, aes(as.factor(mid), y = weight, fill = topo))
  #geom_point(aes(mid, y = weight, fill = topo), shape = 21) +
  #geom_smooth(aes(mid, y = weight, fill = topo)) +
  # geom_segment(aes(x = end, xend = start, y = 1, yend = 1, col = max.topo, group = max.topo), size = 3) +
  # geom_segment(aes(x = start, xend = end, y = 1, yend = 1, col = max.topo, group = max.topo, alpha = as.numeric(topo.max.repeats)), linewidth = 5) +
  geom_segment(aes(x = start, xend = end, y = 0, yend = 0, col = max.topo, group = max.topo, alpha = topo_sig), size = 5) +
  geom_point(aes(mid, y = weight, col = topo), size = 0.2) +
  scale_color_manual(values = topo_cols) +
  facet_grid(chr~., scale = "free_y") +
  theme_classic() 

ggsave(file = paste0(run_name, ".ggplot_all_topos.pdf"), p , width = 10, height = 25)