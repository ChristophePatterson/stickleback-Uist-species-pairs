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
runs <- data.frame(run = gsub(".bestlhoods", "", (list.files()[grepl(analysis_name, list.files())&grepl("\\.bestlhoods", list.files())])))

## Get best module for each pair
runs <- cbind(runs, 
              do.call("rbind", apply(runs, MARGIN = 1, function(x) read.table(paste0(x["run"], ".bestlhoods"), header = T))))

## Write out combined bestlhoods
write.table(runs, paste0(analysis_name, "_all_bootstrap_bestlhoods.txt"), sep = "\t", row.names = F, quote = F)

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

if(analysis_name=="all-monophy-loch_unfolded"){
  r <- ggplot(runs) +
    #### Resident
    # CLAC loch
    geom_segment(aes(x = median(TDivCLAC.), y = "1.0.CLAC", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivCLAC., y = "1.5.")) +
    geom_segment(aes(x = median(TDivCLAC.), y = "2.0.CLAM", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivCLAC.), y = "2.0.CLAM", yend = "1.0.CLAC"), linewidth = 1, lineend = "round") +
    # LUIB loch
    geom_segment(aes(x = median(TDivLUIB.), y = "3.0.LUIB", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivLUIB., y = "3.5.")) +
    geom_segment(aes(x = median(TDivLUIB.), y = "4.0.LUIM", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivLUIB.), y = "4.0.LUIM", yend = "3.0.LUIB"), linewidth = 1, lineend = "round") +
    # DUIN loch
    geom_segment(aes(x = median(TDivDUIN.), y = "5.0.DUIN", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivDUIN., y = "5.5.")) +
    geom_segment(aes(x = median(TDivDUIN.), y = "6.0.DUIM", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivDUIN.), y = "6.0.DUIM", yend = "5.0.DUIN"), linewidth = 1, lineend = "round") +
    # OBSE loch
    geom_segment(aes(x = median(TDivOBSE.), y = "7.0.OBSE", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivOBSE., y = "7.5.")) +
    geom_segment(aes(x = median(TDivOBSE.), y = "8.0.OBSM", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivOBSE.), y = "8.0.OBSM", yend = "7.0.OBSE"), linewidth = 1, lineend = "round") +
    # East lochs
    geom_segment(aes(x = median(TDivOBSE.), y = "7.5.", xend = median(TDivEast.)), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivEast., y = "6.5.")) +
    geom_segment(aes(x = median(TDivDUIN.), y = "5.5.", xend = median(TDivEast.)), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivEast.), y = "7.5.", yend = "5.5."), linewidth = 1, lineend = "round") +
    # Eastlochs
    geom_segment(aes(x = median(TDivLUIB.), y = "3.5.", xend = median(TDivWest.)), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivWest., y = "2.5.")) +
    geom_segment(aes(x = median(TDivCLAC.), y = "1.5.", xend = median(TDivWest.)), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivWest.), y = "3.5.", yend = "1.5."), linewidth = 1, lineend = "round") +
    # MigrAnce
    geom_segment(aes(x = median(TDivAncs.), y = "6.5.", xend = median(TDivEast.)), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivAncs., y = "4.5.")) +
    geom_segment(aes(x = median(TDivAncs.), y = "2.5.", xend = median(TDivWest.)), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivAncs.), y = "6.5.", yend = "2.5."), linewidth = 1, lineend = "round") +
    geom_segment(aes(y = "4.5.", x = median(TDivAncs.), xend = max(TDivAncs.)*1.1), linewidth = 1, lineend = "round")
  
}

if(analysis_name=="all_unfolded"){
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
    geom_segment(aes(x = median(TDivREast.), y = "4.0.OBSE", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivREast.), y = "4.0.OBSE", yend = "3.0.DUIN"), linewidth = 1, lineend = "round") +
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
    geom_segment(aes(x = median(TDivMEast.), y = "8.0.OBSE", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivMEast.), y = "8.0.OBSE", yend = "7.0.DUIM"), linewidth = 1, lineend = "round") +
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
  
}

if(analysis_name=="sigMig"){
  r <- ggplot(runs) +
    #### Resident
    # ResiWest
    geom_segment(aes(x = median(TDivResiWest.), y = "1.0.CLAC", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivResiWest., y = "1.5.")) +
    geom_segment(aes(x = median(TDivResiWest.), y = "2.0.LUIB", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivResiWest.), y = "2.0.LUIB", yend = "1.0.CLAC"), linewidth = 1, lineend = "round") +
    # ResiEast
    geom_segment(aes(x = median(TDivResiEast.), y = "3.0.DUIN", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivResiEast., y = "3.5.")) +
    geom_segment(aes(x = median(TDivResiEast.), y = "4.0.OBSE", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivResiEast.), y = "4.0.OBSE", yend = "3.0.DUIN"), linewidth = 1, lineend = "round") +
    # ResiMerge
    geom_segment(aes(x = median(TDivResi.), y = "1.5.", xend = median(TDivResiWest.)), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivResi., y = "2.5.")) +
    geom_segment(aes(x = median(TDivResi.), y = "3.5.", xend = median(TDivResiEast.)), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivResi.), y = "3.5.", yend = "1.5."), linewidth = 1, lineend = "round") +
    #### Migrate
    # MigrAnce
    geom_segment(aes(x = median(TDivAncs.), y = "6.5.Migr", xend = 0), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivAncs., y = "4.5.")) +
    geom_segment(aes(x = median(TDivAncs.), y = "2.5.", xend = median(TDivResi.)), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivAncs.), y = "6.5.Migr", yend = "2.5."), linewidth = 1, lineend = "round") +
    geom_segment(aes(y = "4.5.", x = median(TDivAncs.), xend = max(TDivAncs.)*1.1), linewidth = 1, lineend = "round")
  
}

if(analysis_name=="sigMandS_unfolded"){
  r <- ggplot(runs) +
    #### Migrate
    # MigrAnce
    geom_segment(aes(x = median(TDivAncs.), y = "1.0.Migr", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivAncs.), y = "2.0.Resi", xend = 0), linewidth = 1, lineend = "round") +
    geom_segment(aes(x = median(TDivAncs.), y = "1.0.Migr", yend = "2.0.Resi"), linewidth = 1, lineend = "round") +
    geom_jitter(aes(x = TDivAncs., y = "1.5.")) +
    geom_segment(aes(y = "1.5.", x = median(TDivAncs.), xend = max(TDivAncs.)*1.1), linewidth = 1, lineend = "round")
  
}


r

# Final tweaks
r <- r + theme_bw() + scale_x_reverse(expand = c(0, 0), name = "Generations (~1 year)")  +
  scale_y_discrete(labels = function(x) substr(x, 5, 8), position = "right",
                   name = "") +
  theme(axis.ticks.y = element_blank())

## Add title to patchwork plot
plot1 <- q / r + plot_annotation(title = paste0(analysis_name))

ggsave(paste0(analysis_name, "bootstrap_combined.png"), plot1, width = 10, height = 10)
