source("~/apps/twisst/plot_twisst.R")
# Read in unique combination of populations
args <- commandArgs(trailingOnly = TRUE)
# set path
top_dir <- args[1]

setwd(top_dir)
source("~/apps/twisst/plot_twisst.R")
library(ggplot2)

pop_combn <- read.table("pop_uniq.txt_combn.txt")
pop_combn$comb <- paste0(pop_combn[,1], "_", pop_combn[,2])
pop_combn <- pop_combn[-grep("OLAV", pop_combn$comb),]

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
scaf <- read.table("GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
twisst_data_all$chr <- scaf$Sequence.name[match(twisst_data_all$scaffold,scaf$RefSeq.seq.accession)]

p <- ggplot(twisst_data_all) +
  geom_step(aes(start, topo2, col = run_name), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x") +
  # scale_color_manual(values = topo_cols) +
  theme_classic()

ggsave(file = paste0(top_dir, "/All_population_comparison_twisst.pdf"), p, width = 15, height = 30)
ggsave(file = paste0(top_dir, "/All_population_comparison_twisst.png"), p, width = 15, height = 30)


