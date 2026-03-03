library(tidyverse)

## Get commnad arguments
args <- commandArgs(trailingOnly = TRUE)

# Get data
inputdir <- args[1]

# Read in table
cov_df <- read.table(paste0(inputdir, "/Sex_coverage.txt"), header = T) %>%
  tibble()

# Divide X coverage by average autosomal coverage for each individual
cov_df <- cov_df %>%
  group_by(individual) %>%
  mutate(cov.dv.mean = meandepth/mean(meandepth[rname=="CM102076.1"])) %>%
  ungroup()

p <- ggplot(cov_df) +
  geom_vline(xintercept = 2500000) +
  geom_vline(xintercept = 12500000) +
  geom_line(aes(startpos, cov.dv.mean, col = individual), show.legend = F) +
  facet_grid(rname~., space = "free") +
  theme_bw() +
  coord_cartesian(ylim = c(0,2)) +
  scale_x_continuous(n.breaks = 10, labels = function(x) paste0(x / 1e6))
ggsave(paste0(inputdir, "/Sex_coverage_across_genome.png"), p, width = 15, height = 5)

# Calculate average X and Y coverage
cov_sum <- cov_df %>% 
  group_by(individual, .add = T) %>%
  summarise(mn.auto.depth = mean(meandepth[rname=="CM102076.1"], na.rm = T),
            mn.PAR.depth = mean(meandepth[rname=="CM102094.1"&between(startpos,0, 2500000)], na.rm = T),
            mn.XY.depth = mean(meandepth[rname=="CM102094.1"&between(startpos,2500000,12500000)], na.rm = T),
            mn.X.depth = mean(meandepth[rname=="CM102094.1"&startpos>12500000], na.rm = T))

## Set ratio limit for sex determination
Xcov_sex_determine <- 0.75

# Calculate sex
cov_sum$sex.ratio <- cov_sum$mn.X.depth/cov_sum$mn.auto.depth
cov_sum$Gsex <- NA
cov_sum$Gsex[cov_sum$sex.ratio>=Xcov_sex_determine] <- "F"
cov_sum$Gsex[cov_sum$sex.ratio<=Xcov_sex_determine] <- "M"

# Create bed file
bedfile <- data.frame(fam = "Fam1", individual = cov_sum$individual, X = 0, Y = 0, Gsex = cov_sum$Gsex)

# Write out bed file
write.table(bedfile, quote = F,row.names = F,col.names = F, file = paste0(inputdir, "/Gsex.ped"), sep = "\t")

## Plot ratio of X to Y
p <- ggplot(cov_sum) +
  geom_jitter(aes(Gsex, mn.X.depth/mn.auto.depth, col = Gsex), size = 3, height = 0) +
  geom_jitter(aes(Gsex, mn.X.depth/mn.auto.depth, col = Gsex), size = 3, height = 0) +
  lims(y = c(0, NA)) +
  theme_bw()

ggsave(paste0(inputdir, "/Sex_coverage.png"), p)