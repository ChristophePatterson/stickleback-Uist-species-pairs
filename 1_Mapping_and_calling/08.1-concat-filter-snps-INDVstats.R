library(tidyverse)
library(patchwork)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf <- args[1]
wkdir <- dirname(vcf)

stat_ext <- paste0(wkdir, "/stats/", gsub(".bcf", "", basename(vcf)))
stat_ext <- gsub(".vcf.gz", "", stat_ext)

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

## Get estimates of sex from depth ofX chromosome
Gsex <- read_table("/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Gsex.ped", col_names = F) %>% 
rename(individual= X2,gsex = X5) %>%
select(individual, gsex)

idf <- tibble(merge(idf, Gsex, by = "individual"))

# Calculate prop of heterozygous samples
idf$prop_het <- 1-(idf$OHOM/idf$N_SITES)

## Calculate excess heterozygousity
idf <- idf %>%
  group_by(Population) %>%
  mutate(excess_het = prop_het-(mean(prop_het)),
         excess_het_sd = ((prop_het-(mean(prop_het)))>0.015)) %>%
  ungroup()

# Plot
p <- ggplot(idf) +
  geom_boxplot(aes(Population, MEAN_DEPTH), outlier.shape = NA) +
  geom_jitter(aes(Population, MEAN_DEPTH, col = Population), height = 0, width = 0.2) +
  lims(y = c(0, 50)) +
  ggtitle("Mean sample depth vs pop") +
  scale_y_continuous(name = "Mean Depth") + theme_bw()

q <- ggplot(idf) +
  geom_boxplot(aes(Population, N_MISS/N_DATA), outlier.shape = NA) +
  geom_jitter(aes(Population, N_MISS/N_DATA, col = Population), height = 0, width = 0.2) +
  ggtitle("Sample missing data vs pop") +
  scale_y_continuous(name = "Proportion of Missing Sites") + theme_bw()

r <- ggplot(idf) +
  geom_point(aes(as.numeric(sum_file_size), N_MISS/N_DATA, col = Population)) +
  scale_x_log10() +
  ggtitle("Sample seq-data size vs missing data") +
  scale_y_continuous(name = "Proportion of Missing Sites") +
  scale_x_continuous(name = "Total sequence data (log10)") + theme_bw() 

s <- ggplot(idf) +
  geom_point(aes(N_MISS/N_DATA, OHOM/N_SITES, col = Population), height = 0, width = 0.2) +
  ggtitle("Sample homozgousity vs missing data") + 
  scale_y_continuous(name = "Proportion of homozygous sites") +
  scale_x_continuous(name = "Number of missing sites") +
  theme_bw()

s1 <- ggplot(idf) +
  geom_boxplot(aes(Population, prop_het), outlier.shape = NA) +
  geom_jitter(aes(Population, prop_het), height = 0, width = 0.2, size = 3) +
  ggtitle("Sample homozgousity by population") +
  scale_y_continuous(name = "Proportion of heterozygous sites") + theme_bw()

## Save plot
plot1 <- p+q+r+s + plot_layout(guides = "collect")

paste0(stat_ext,"stats_plot.png")
ggsave(paste0(stat_ext,"stats_plot.png"), plot1, width = 10, height = 10)
ggsave(paste0(stat_ext,"stats_het_pop_plot.png"), s1, width = 8, height = 8)

## Get samples that are high QUAL
HQsamples <- idf$individual[idf$N_MISS/idf$N_DATA<=0.4&idf$excess_het<0.015]
# Also remove sample NOVSC116 which contains some odd mitochondrial data
HQsamples <- HQsamples[HQsamples!="NOVSC116"]

print(paste0(sum(idf$N_MISS/idf$N_DATA>=0.4), " samples do not meet critearia for missing data"))
print(paste0(sum(idf$excess_het>0.015), " samples have excess heterozgousity"))
writeLines(con = paste0(stat_ext, "_HiQ_vcf_samples.txt"), HQsamples)

## Save merged stats file
idf
write.table(idf, paste0(stat_ext,"_merged_stats.txt"), row.names = F, sep = "\t")

