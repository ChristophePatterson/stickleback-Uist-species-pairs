library(tidyverse)
library(patchwork)
library(ggnewscale)

# Get command line arguments
args <- commandArgs(trailingOnly = T)

# Set working directory
setwd(args[1])

# Get analysis name
analysis_name <- args[2]

# Get list of pop pairs
runs <- data.frame(run = gsub(".bestlhoods", "", list.files()[grep(".bestlhoods", list.files())]))

# Get population names
pop.order <- c("CLAC","CLAM", "DUIN", "DUIM", "LUIB","LUIM", "OBSE", "OBSM")

## Get best module for each pair
runs <- cbind(runs, 
                   do.call("rbind", apply(runs, MARGIN = 1, function(x) read.table(paste0(x["run"], ".bestlhoods"), header = T))))

runs_long <- pivot_longer(runs, cols = colnames(runs)[-1]) %>%
  mutate(name = factor(name, levels = colnames(runs)),
         is.div = grepl("TDiv", name),
         is.popNe = !grepl("TDiv|Lhood",name))
  

## Plot of best pop sizes
p <- ggplot(runs_long[runs_long$is.div,],aes(name, value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2)

q <- ggplot(runs_long[runs_long$is.popNe,],aes(name, value/2)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  scale_x_discrete(labels = function(x) gsub("\\.","", x), name = "Population") +
  scale_y_continuous(name = "Ne (1N in Thousands)", labels = function(x) x/1000) +
  theme_bw()

distri_plot <- p + q & theme_bw()

r <- ggplot(runs) +
  #### Resident
  # ResiWest
  geom_segment(aes(x = median(TDivRWest.), y = "1.0.CLAC", xend = 0), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivRWest., y = "1.5.")) +
  geom_segment(aes(x = median(TDivRWest.), y = "2.0.LUIB", xend = 0), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivRWest.), y = "2.0.LUIB", yend = "1.0.CLAC"), linewidth = 1, lineend = "round") +
  # ResiEast
  geom_segment(aes(x = median(TDivREast.), y = "3.0.DUIN", xend = 0), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivREast., y = "3.5.")) +
  geom_segment(aes(x = median(TDivREast.), y = "4.0.LUIB", xend = 0), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivREast.), y = "4.0.LUIB", yend = "3.0.DUIN"), linewidth = 1, lineend = "round") +
  # ResiMerge
  geom_segment(aes(x = median(TDivResi.), y = "1.5.", xend = median(TDivRWest.)), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivResi., y = "2.5.")) +
  geom_segment(aes(x = median(TDivResi.), y = "3.5.", xend = median(TDivREast.)), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivResi.), y = "3.5.", yend = "1.5."), linewidth = 1, lineend = "round") +
  #### Migrate
  # MigrWest
  geom_segment(aes(x = median(TDivMWest.), y = "5.0.CLAM", xend = 0), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivMWest., y = "5.5.")) +
  geom_segment(aes(x = median(TDivMWest.), y = "6.0.LUIM", xend = 0), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivMWest.), y = "6.0.LUIM", yend = "5.0.CLAM"), linewidth = 1, lineend = "round") +
  # MigrEast
  geom_segment(aes(x = median(TDivMEast.), y = "7.0.DUIM", xend = 0), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivMWest., y = "7.5.")) +
  geom_segment(aes(x = median(TDivMEast.), y = "8.0.LUIM", xend = 0), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivMEast.), y = "8.0.LUIM", yend = "7.0.DUIM"), linewidth = 1, lineend = "round") +
  # MigrMerge
  geom_segment(aes(x = median(TDivMigr.), y = "7.5.", xend = median(TDivMWest.)), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivMigr., y = "6.5.")) +
  geom_segment(aes(x = median(TDivMigr.), y = "5.5.", xend = median(TDivMEast.)), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivMigr.), y = "7.5.", yend = "5.5."), linewidth = 1, lineend = "round") +
  # MigrAnce
  geom_segment(aes(x = median(TDivAncs.), y = "6.5.", xend = median(TDivMigr.)), linewidth = 1, lineend = "round") +
  geom_jitter(aes(x = TDivAncs., y = "4.5.")) +
  geom_segment(aes(x = median(TDivAncs.), y = "2.5.", xend = median(TDivResi.)), linewidth = 1, lineend = "round") +
  geom_segment(aes(x = median(TDivAncs.), y = "6.5.", yend = "2.5."), linewidth = 1, lineend = "round") +
  geom_segment(aes(y = "4.5.", x = median(TDivAncs.), xend = max(TDivAncs.)*1.1), linewidth = 1, lineend = "round")

# Final tweaks
r <- r + theme_bw() + scale_x_reverse(expand = c(0, 0), name = "Generations (~1 year)")  +
  scale_y_discrete(labels = function(x) substr(x, 5, 8), position = "right",
                  name = "") +
  theme(axis.ticks.y = element_blank())

ggsave(paste0(analysis_name, "bootstrap_combined.png"), q / r, width = 10, height = 10)
