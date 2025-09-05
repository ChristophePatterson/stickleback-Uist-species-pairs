library(tidyverse)

ROH <- read_table("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/All.RG.txt", skip = 1, 
                    col_names = c("RG","Sample","Chromosome","Start","End","Length","Number.of.markers", "Quality"))

G.length <- as.numeric(readLines("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/GCA_046562415.1_Duke_GAcu_1.0_genomic.length.txt"))

ROH.sum <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/het_counts.txt", header = T)

ROH.sum <- ROH.sum %>%
  mutate(nonROHhet = hetn/(G.length-ROHsum)*1000,
         FROH = ROHsum/G.length) %>%
  mutate(IBrisk = nonROHhet*FROH)

ROH.pop <- ROH.sum %>%
  group_by(pop) %>%
  summarise(mn.FROH = mean(FROH),
    mn.nonROHhet = mean(nonROHhet),
    IBrisk = mn.FROH*mn.nonROHhet)

p <- ggplot(ROH.pop)+
 geom_point(data = ROH.sum, aes(FROH, nonROHhet, color = pop)) +
 geom_text(aes(mn.FROH, mn.nonROHhet, label = pop)) 

ggsave("test.png", p)


p <- ggplot(ROH)+
 geom_histogram(aes(Length)) +
 facet_wrap(~Sample)

ggsave("test.png", p)

## Get range of cutoff lengths for ROH
params <- expand_grid(
  group = unique(ROH$Sample),
  cutoff = seq(0, 20e5, 1e5)
)

# Calculate Frequence of ROH for each sample across the genome
FROH <- params %>%
  mutate(stats = map2(group, cutoff, ~ {
    ROH %>%
      filter(Sample == .x, Length > .y) %>%
      summarise(mean_val = mean(Length, na.rm = TRUE),
                FROH = sum(Length)/G.length,
                n        = n())
  })) %>%
  unnest(stats) %>%
  rename(Sample = group)

## Read in sample data
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
FROH <- merge(FROH, samples, by.x = "Sample", by.y="individual", all.x = T)

p <- ggplot(FROH) +
    geom_boxplot(aes(Population, FROH)) +
    facet_wrap(~cutoff)

ggsave("test.png", p)


p.ROH <- ggplot(ROH) +
    geom_segment(aes(x = Start, xend = End, y = Sample, yend = Sample, color = Sample), show.legend = F) +
    facet_grid(Chromosome~.)

ggsave("test.png", p.ROH , width = 10, height = 20)