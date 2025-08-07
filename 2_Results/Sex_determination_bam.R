library(tidyverse)

## Get commnad arguments
args <- commandArgs(trailingOnly = TRUE)

# Get data
inputdir <- args[1]

# Read in table
cov_df <- read.table(paste0(inputdir, "/Sex_coverage.txt"), header = T)

# Calculate average X and Y coverage
cov_sum <- cov_df %>% 
  group_by(individual, .add = T) %>%
  summarise(mn.auto.depth = mean(meandepth[rname=="CM102076.1"], na.rm = T),
            mn.X.depth = mean(meandepth[rname=="CM102094.1"], na.rm = T))

## Set ratio limit for sex determination
Xcov_sex_determine <- 0.85

# Calculate sex
cov_sum$sex.ratio <- cov_sum$mn.X.depth/cov_sum$mn.auto.depth
cov_sum$Gsex <- NA
cov_sum$Gsex[cov_sum$sex.ratio>=Xcov_sex_determine] <- "F"
cov_sum$Gsex[cov_sum$sex.ratio<=Xcov_sex_determine] <- "M"

# Create bed file
bedfile <- data.frame(individual = cov_sum$individual, X = 0, Y = 0, Gsex = cov_sum$Gsex)

# Write out bed file
write.table(bedfile, quote = F,row.names = F,col.names = F, file = paste0(inputdir, "/Gsex.bed"), sep = "\t")

## Plot ratio of X to Y
p <- ggplot(cov_sum) +
  geom_point(aes(individual, mn.X.depth/mn.auto.depth, col = Gsex), size = 3)

ggsave(paste0(inputdir, "/Sex_coverage.png"), p)