library(tidyverse)
library(patchwork)

## Colorblind palette
cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
ROH.file <- args[1]
# ROH.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/All.RG.txt"
ROH.calc.file <- args[2]
# ROH.calc.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/All.RG.calcs.txt"

# Read in position of all runs of heterozgousity
ROH <- read_table(ROH.file, skip = 1, 
                    col_names = c("RG","Sample","Chromosome","Start","End","Length","Number.of.markers", "Quality"))

# Get total length of genome
G.length <- as.numeric(readLines("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/GCA_046562415.1_Duke_GAcu_1.0_genomic.length.txt"))

## Get ROH stats calcs
ROH.calcs <- read.csv(ROH.calc.file, header = T)

## Read in sample data
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
ROH.calcs <- merge(ROH.calcs, samples, by.x = "samp", by.y="individual", all.x = T)

## Calculations
ROH.calcs <- ROH.calcs %>%
  mutate(nonROHhet = nHet/(G.length-sumROH)*1000, # nonROH heterozgousity as per 1kb
         FROH = sumROH/G.length) # Proportion of genome that is identified as ROH

# Calculate popuation statistics
ROH.pop <- ROH.calcs %>%
  group_by(Population) %>% # group by population
  group_by(cutoff, .add = T) %>% # group by ROH cutoff
  summarise(mn.FROH = mean(FROH), # Mean of FROH
    mn.nonROHhet = mean(nonROHhet), # Mean of nonROH heterozygousity 
    IBrisk = mn.FROH*mn.nonROHhet) # IBrisk per population

## PLot indv samples across cutoff
p <- ggplot(ROH.calcs)+
  geom_line( aes(cutoff, FROH, color = Population, group = samp)) +
  geom_point(aes(cutoff, FROH, color = Population)) +
  geom_text(aes(cutoff, FROH, color = Population, label = nROH), vjust = -0.5, show.legend = F) +
  ylab(bquote(F[ROH])) +
  scale_color_manual(values = cbPalette) +
  # facet_grid(Population~.) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6)) +
  xlab("ROH cutoff (Mbps)") +
  theme_bw()


ggsave("test.png", p, width = 10)

q <- ggplot(ROH.pop)+
  geom_line(aes(cutoff, IBrisk, color = Population, group = Population)) +
  geom_point(aes(cutoff, IBrisk, fill = Population), shape = 21, size = 2) +
  ylab(bquote(IB[risk])) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  # facet_grid(Population~.) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6)) +
  xlab("ROH cutoff (Mbps)") +
  theme_bw()

ggsave(paste0(gsub(".txt", "", ROH.calc.file),".FROH_IBrisk.png"), p/q, width = 10)

### Plot position of ROHs
ROH <- merge(ROH, samples[,c("individual", "Population", "Ecotype")], by.x = "Sample", by.y="individual", all.x = T)

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
ROH$chr <- factor(as.character(as.roman(chr$Chromosome.name[match(ROH$Chromosome, chr$GenBank.seq.accession)])), levels = as.character(as.roman(1:22)))


p.ROH <- ggplot(ROH) +
    geom_segment(aes(x = Start, xend = End, y = Sample, yend = Sample, color = Population), show.legend = F, linewidth = 3) +
    facet_grid(Population~chr, space = "free", scales = "free") +
    scale_color_manual(values = cbPalette) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6)) +
    theme_bw() +
    theme(panel.spacing = unit(0, "points", data = NULL))

ggsave(paste0(gsub(".txt", "", ROH.calc.file),".ROHs.png"), p.ROH, width = 20, height = 10)