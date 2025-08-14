library(tidyverse)
library(patchwork)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf <- args[1]
wkdir <- dirname(vcf)

stat_ext <- paste0(wkdir, "/stats/", gsub(".vcf.gz", "", basename(vcf)))

# Directory location
setwd(wkdir)

# Paired waterbodies
paired_waterbodies <- c("CLAC", "OBSE", "LUIB", "DUIN")

# Read in het
ihet <- read_table(paste0(stat_ext, ".het")) %>%
  rename_with(~ gsub("\\(|\\)", "", .x))
# Read in missing data
imiss <- read_table(paste0(stat_ext, ".imiss")) %>%
  rename_with(~ gsub("\\(|\\)", "", .x))

# Read in depth
idepth <- read_table(paste0(stat_ext, ".idepth")) %>%
  rename_with(~ gsub("\\(|\\)", "", .x))

# Merge all into one
idf <- tibble(merge(merge(ihet, imiss , by = "INDV"), idepth[,colnames(idepth)!="N_SITES"]))

### Combine Uist samples with Uist sample locations
seq_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
idf <- tibble(merge(seq_data, idf, by.x = "individual", by.y = "INDV"))

# Plot
p <- ggplot(idf) +
  geom_boxplot(aes(Population, MEAN_DEPTH), outlier.shape = NA) +
  geom_jitter(aes(Population, MEAN_DEPTH, col = Population), height = 0, width = 0.2) +
  lims(y = c(0, 50)) +
  ggtitle("Mean sample depth vs pop")

q <- ggplot(idf) +
  geom_boxplot(aes(Population, N_MISS/N_DATA), outlier.shape = NA) +
  geom_jitter(aes(Population, N_MISS/N_DATA, col = Population), height = 0, width = 0.2) +
  ggtitle("Sample missing data vs pop")

r <- ggplot(idf[idf$Waterbody%in%paired_waterbodies,]) +
  geom_point(aes(as.numeric(sum_file_size), N_MISS/N_DATA, col = Population)) +
  scale_x_log10() +
  ggtitle("Sample seq-data size vs missing data")

s <- ggplot(idf[idf$Waterbody%in%paired_waterbodies,]) +
  geom_point(aes(N_MISS, OHOM/N_SITES, col = Population), height = 0, width = 0.2) +
  ggtitle("Sample homozgousity vs missing data")

## Save plot
plot1 <- p+q+r+s + plot_layout(guides = "collect")

ggsave("stats/stats_plot.png", plot1, width = 16, height = 12)

## Get samples that are high QUAL
HQsamples <- idf$individual[idf$N_MISS/idf$N_DATA<=0.4]

print(paste0(sum(idf$N_MISS/idf$N_DATA>=0.4), " samples do not meet critearia for missing data"))
writeLines(con = "HiQ_vcf_samples.txt", HQsamples)

