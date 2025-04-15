source("~/apps/twisst/plot_twisst.R")

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


library(tidyverse)

# Convert twisst data into tidy format
window_data_df <- do.call("rbind", twisst_data[["window_data"]])
weight_data_df <- do.call("rbind", twisst_data[["weights_raw"]])
# Combine window and weight data
twisst_data_sf <- cbind(window_data_df, weight_data_df)
# Convert into tidy format
twisst_data_sf <- pivot_longer(twisst_data_sf, cols = paste0("topo", 1:3), names_to = "topo", values_to = "weight")
window_ave <- mean(twisst_data_sf$end - twisst_data_sf$start, na.rm = T)

### 
pdf(file = paste0(run_name, ".summary.pdf"))
plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
dev.off()

pdf(file = paste0(run_name, ".all_typo.pdf"), width  = 10, height = 40)
plot.twisst(twisst_data, mode=1, show_topos=TRUE)
dev.off()

pdf(file = paste0(run_name, ".all_typo_stacked.pdf"), width  = 10, height = 40)
plot.twisst(twisst_data, mode=3, show_topos=TRUE)
dev.off()
### 
### # make smooth weightings and plot those across chromosomes
chr <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
chr_select <- chr$RefSeq.seq.accession[chr$Seq.length>10000000&!grepl("Y", chr$Chromosome.name)]
twisst_data_subset <- subset.twisst.by.regions(twisst_data, regions = chr_select)


# Replace scaffold name with chr name
twisst_data_sf$chr <- chr$Chromosome.name[match(twisst_data_sf$scaffold, chr$RefSeq.seq.accession)]


p <- ggplot(twisst_data_sf) +
  #geom_bar(stat="identity", position = "stack", width = 1, aes(as.factor(mid), y = weight, fill = topo))
  #geom_point(aes(mid, y = weight, fill = topo), shape = 21) +
  #geom_smooth(aes(mid, y = weight, fill = topo)) +
  geom_step(aes(start, y = weight, col = topo)) +
  facet_grid(chr~.) +
  theme_classic()

ggsave(file = paste0(run_name, ".ggplot_all_topos.pdf"), p , width = 10, height = 25)