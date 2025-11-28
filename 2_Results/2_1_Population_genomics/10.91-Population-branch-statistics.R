library(tidyverse)
library(patchwork)
library(tidyverse)
library(ggnewscale)

args <- commandArgs(trailingOnly = TRUE)

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# set path
my_bins <- args[1]

## Set order pops should be plotted
allpops <- c("resi", "anad", "fw")

# Read in data
GenomicPops <- read_csv(paste0(my_bins,"_APARX.csv"))
GenomicPops <- read_csv("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/PopBranchStats/sliding_window_w25kb_s5kb_m1_PopPair_APARX.csv")

## Convert to distance along genome
scaf <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

# Add chr name to GenomicPops
GenomicPops$chr <- scaf$Sequence.name[match(GenomicPops$scaffold, scaf$GenBank.seq.accession)]
# Calculate cumulative length (starting at 0 )
scaf$Cum.Seq.length <- c(0, cumsum(scaf$Seq.length[1:(nrow(scaf)-1)]))
# Add cumulative length to GenomicPops
GenomicPops$start.cum <- GenomicPops$start+(scaf$Cum.Seq.length[match(GenomicPops$chr, scaf$Sequence.name)]+scaf$Cum.Seq.length[1])



## Convert fst into -log(1-Fst)
GenomicPops <- GenomicPops %>%
  mutate(across(starts_with("Fst_"), ~ -log10(1 - .), .names = "T_{.col}"))

# Calculate PBS for each population
GenomicPops$PBS_resi <- (GenomicPops$T_Fst_resi_fw + GenomicPops$T_Fst_resi_anad - GenomicPops$T_Fst_anad_fw) / 2
GenomicPops$PBS_fw <- (GenomicPops$T_Fst_resi_fw + GenomicPops$T_Fst_anad_fw - GenomicPops$T_Fst_resi_anad) / 2
GenomicPops$PBS_anad <- (GenomicPops$T_Fst_resi_anad + GenomicPops$T_Fst_anad_fw - GenomicPops$T_Fst_resi_fw) / 2
# Reshape data for plotting
GenomicPops_long <- GenomicPops %>%
  select(scaffold, chr, start, end, start.cum, starts_with("PBS_")) %>%
  pivot_longer(cols = starts_with("PBS_"), names_to = "Population", values_to = "PBS") %>%
  mutate(Population = str_replace(Population, "PBS_", ""))

# regions <- data.frame(chr = factor(c("I", "IV", "XI", "XXI"), levels = levels(chr$Sequence.name)), start = c(25000000, 12000000, 5000000, 8000000), end = c(31000000, 16000000, 10000000, 15000000))
regions <- data.frame(chr = factor(c("chrI", "chrIV", "chrIX", "chrXI", "chrXIX", "chrXXI")), start = c(26000000, 12000000, 12000000, 5500000, 2000000, 8000000), end = c(27500000, 16000000, 14500000,7000000, 10000000, 15000000))

# Filer dataset to specific region
GenomicPops_long_filt <- GenomicPops_long %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()


# Rename populations
GenomicPops_long$Population[GenomicPops_long$Population=="anad"] <- "mig"

# Plot PBS
pbs_plot <- ggplot(GenomicPops_long, aes(x = start.cum, y = PBS, color = Population)) +
  geom_vline(xintercept = scaf$Cum.Seq.length) +
  geom_line() +
  # Add vertical lines for chromosome breaks
  # Add chromosome labels
  scale_x_continuous(breaks = (scaf$Cum.Seq.length + scaf$Seq.length / 2),
                     labels = scaf$Sequence.name) +
  scale_color_manual(values = cbPalette[1:length(allpops)]) +
  labs(x = "Genomic Position", y = "Population Branch Statistic (PBS)") +
  theme_bw() +
  theme(legend.position = "bottom")    +
  facet_wrap(~ Population, ncol = 1) +
  theme(strip.text = element_text(size = 12, face = "bold")) 

# For filtered regions
pbs_regions <- ggplot(GenomicPops_long_filt, aes(x = start, y = PBS, color = Population)) +
  geom_line() +
  scale_color_manual(values = cbPalette[1:length(allpops)]) +
  labs(x = "Genomic Position", y = "Population Branch Statistic (PBS)") +
  theme_bw() +
  theme(legend.position = "bottom")    +
  facet_grid(Population ~ chr, scales = "free_x") +
  theme(strip.text = element_text(size = 12, face = "bold")) 

ggsave("test.png", plot = pbs_plot/pbs_regions, width = 20, height = 15)

## Population Fst long format
FstPops_long <- pivot_longer(GenomicPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("Fst_", colnames(GenomicPops))])], 
                             cols = colnames(GenomicPops)[grep("Fst_", colnames(hetPops))], values_to = "Fst", names_prefix = "Fst_", names_to = "Population")


fst_diff <- ggplot(GenomicPops, aes(x = start, y = T_Fst_resi_fw + T_Fst_resi_anad - T_Fst_anad_fw)) +
    geom_line() +
    geom_line(aes(y = -Fst_anad_fw, col = "anad_fw")) +
    geom_line(aes(y = +Fst_resi_anad, col = "resi_fw")) +
    labs(x = "Genomic Position", y = "Fst(resi-mig) - Fst(mig-fw)") +
      scale_color_manual(values = cbPalette) +
    facet_grid(chr~.) +
    theme_bw() +
    theme(legend.position = "bottom")


ggsave("test.png", plot = fst_diff, width = 20, height = 15)


fst_diff <- ggplot(GenomicPops, aes(Fst_anad_fw, Fst_resi_anad)) +
    geom_point() +
    facet_wrap(~chr) +
    theme_bw() +
    theme(legend.position = "bottom")


ggsave("test.png", plot = fst_diff, width = 20, height = 15)