# Output file location
library(tidyverse)
library(patchwork)
args <- commandArgs(trailingOnly = TRUE)

# set path
SNP.library.name <- args[1]
my_bins <- args[2]
print(my_bins)

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
p <- ggplot(ld_bins, aes(distance, avg_R2)) + geom_point(alpha = 0.3) + geom_smooth() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))

q <- p + xlim(c(0, 100000))

# Save file
ggsave(paste0(SNP.library.name, ".LD.ld_decay_bins.pdf"), q + p, height = 10, width = 20)
# ggsave(paste0("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/", SNP.library.name, ".LD.ld_decay_bins.pdf"), p, height = 10, width = 20)
