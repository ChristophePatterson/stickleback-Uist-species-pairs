library(tidyverse)
library(patchwork)
library(tidyverse)
library(ggnewscale)

args <- commandArgs(trailingOnly = TRUE)

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# set path
my_bins <- args[1]
# my_bins <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/indPops/sliding_window_w100kb_s100kb_m1_PopPair_auto"

## Set order pops should be plotted
allpops <- c("CLAC" ,"DUIN","LUIB","OBSE",  "CLAM", "DUIM", "LUIM", "OBSM")

sex_chrom <- "CM102094.1"

# Read in data
hetPops <- read_csv(paste0(my_bins,"_APARX.csv"))

## Pivot pi
hetPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("pi_", colnames(hetPops))])], 
             cols = colnames(hetPops)[grep("pi_", colnames(hetPops))], values_to = "pi", names_prefix = "pi_", names_to = "Population")
## Pivot dxy
dxyPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("dxy_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("dxy_", colnames(hetPops))], values_to = "dxy", names_prefix = "dxy_", names_to = "Population")
## Pivot fst
FstPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("Fst_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("Fst_", colnames(hetPops))], values_to = "Fst", names_prefix = "Fst_", names_to = "Population")
## Pivot hetInd
# Read in data
hetInd <- read_csv(paste0(my_bins,"_auto.csv"))
HetInd_long <- pivot_longer(hetInd[, c("scaffold","start","end","mid", "sites",colnames(hetInd)[grep("het_", colnames(hetInd))])], 
                             cols = colnames(hetInd)[grep("het_", colnames(hetInd))], values_to = "het", names_prefix = "het_", names_to = "individual")

## Summarise dxy
dxyPops_long_summary <- dxyPops_long %>%
  filter(scaffold != sex_chrom) %>%
  group_by(by = Population) %>%
  mutate(Pop1 = str_split_i(Population, "_", 1),
         Pop2 = str_split_i(Population, "_", 2)) %>%
  summarise(dxy.mn = mean(dxy, na.rm =T), dxy.sd = sd(dxy), 
            Pop1 = Pop1[1], Pop2 = Pop2[1])

## Order pop1 and pop2 so they are the same order as that in allpops
dxyPops_long_summary <- apply(dxyPops_long_summary, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops))]
}
      ) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(dxyPops_long_summary[,1:3])

## Summarise Fst
FstPops_long_summary <- FstPops_long %>%
  filter(scaffold != sex_chrom) %>%
  group_by(by = Population) %>%
  mutate(Pop1 = str_split_i(Population, "_", 1),
         Pop2 = str_split_i(Population, "_", 2)) %>%
  summarise(Fst.mn = mean(Fst, na.rm =T), Fst.sd = sd(Fst, na.rm = T), 
            Pop1 = Pop1[1], Pop2 = Pop2[1])

## Order pop1 and pop2 so they are the same order as that in allpops
FstPops_long_summary <- apply(FstPops_long_summary, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops), decreasing = T)]
}
) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(FstPops_long_summary[,1:3])

## Set factor level and reset POp1 and Pop2
dxyPops_long_summary$Pop1 <- factor(dxyPops_long_summary$V1,levels =  allpops)
dxyPops_long_summary$Pop2 <- factor(dxyPops_long_summary$V2,levels =  allpops)

## Set factor level and reset POp1 and Pop2
FstPops_long_summary$Pop1 <- factor(FstPops_long_summary$V1,levels =  allpops)
FstPops_long_summary$Pop2 <- factor(FstPops_long_summary$V2,levels =  allpops)

# Get summary of pi
hetPops_long_summary <- hetPops_long %>%
  filter(scaffold != sex_chrom) %>%
  group_by(by = Population) %>%
  summarise(pi.mn = mean(pi, na.rm =T), pi.sd = sd(pi, na.rm = T), pi.md = median(pi, na.rm =T)) %>%
  rename(Pop = by)

hetPops_long_summary$Pop <- factor(hetPops_long_summary$Pop, levels = allpops)

pi.plot <- ggplot(hetPops_long) +
  geom_boxplot(aes(Population, pi), outlier.shape = NA) +
  theme_classic()
  
ggsave(paste0(my_bins, "_pi_by_pop.png"), pi.plot, height = 10, width = 10.5)
ggsave(paste0(my_bins, "_pi_by_pop.pdf"), pi.plot, height = 10, width = 10.5)


# Write out stats table
write.table(x = hetPops_long_summary, file = paste0(my_bins, "_hetPops_stats.txt"))

## plot
p <- ggplot() +
  geom_tile(data = dxyPops_long_summary, aes(Pop1, Pop2, fill = dxy.mn)) +
  scale_fill_gradient(low = "white", high = "deepskyblue",name = "Dxy") +
  new_scale_fill() +
  geom_tile(data = FstPops_long_summary, aes(Pop1, Pop2, fill = Fst.mn)) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "orange",name = "Fst") +
  new_scale_fill() +
  geom_tile(data = hetPops_long_summary, aes(Pop, Pop, fill = pi.mn)) +
  scale_fill_gradient(low = "white", high = "firebrick",name = "pi") +
  geom_vline(xintercept = 4.5) +
  geom_hline(yintercept = 4.5) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "grey")) +
  ggtitle("Dxy, pi, and Fst - Lagoon species pairs")

ggsave(paste0(my_bins, ".png"), p, height = 10, width = 10.5)
ggsave(paste0(my_bins, ".pdf"), p, height = 10, width = 10.5)

### Plot of Heterozgousity
## Add in sample data
sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
HetInd_long$Population <- sample_data$Population[match(HetInd_long$individual, sample_data$individual)]
HetInd_long$Population <- factor(HetInd_long$Population, levels = allpops)
HetInd_long$Ecotype <- sample_data$Ecotype[match(HetInd_long$individual, sample_data$individual)]
HetInd_long$Waterbody <- sample_data$Waterbody[match(HetInd_long$individual, sample_data$individual)]

##
HetInd.summary <- HetInd_long %>%
  group_by(by=individual) %>%
  summarise(mn.het = mean(het, na.rm = T), Population = first(Population), Ecotype = first(Ecotype))

phet <- ggplot(HetInd.summary, aes(Population, mn.het)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.1)

ggsave(paste0(my_bins, "_ave_het.png"), phet)


## Plot Fst across the genome
# Check if Fst and Pdxy can be combined
any(!FstPops_long$scaffold==dxyPops_long$scaffold)
any(!FstPops_long$start==dxyPops_long$start)

# Combine
Pops_long <- cbind(FstPops_long, dxyPops_long[,"dxy"])
# Extract first and second population
Pops_long <- Pops_long %>%
  mutate(pop1 = str_split_i(string = Population, pattern = "_", i = 1 ),
         pop2 = str_split_i(string = Population, pattern = "_", i = 2 ))

# Read in chr info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

# Replace scaffold name with chr name
Pops_long$chr <- gsub("chr","",chr$Sequence.name[match(Pops_long$scaffold, chr$GenBank.seq.accession)])
Pops_long$chr <- factor(Pops_long$chr, levels = gsub("chr", "", chr$Sequence.name[order(chr$GenBank.seq.accession)]))

# Determine Ecotype comparison
Pops_long$Ecotypes <- NA
Pops_long$Ecotypes[Pops_long$pop1%in%(allpops[1:4])&Pops_long$pop2%in%(allpops[1:4])] <- "resi_resi"
Pops_long$Ecotypes[Pops_long$pop1%in%(allpops[5:8])&Pops_long$pop2%in%(allpops[5:8])] <- "mig_mig"
Pops_long$Ecotypes[(Pops_long$pop1%in%(allpops[1:4])&Pops_long$pop2%in%(allpops[5:8]))] <- "resi_mig"
Pops_long$Ecotypes[(Pops_long$pop1%in%(allpops[5:8])&Pops_long$pop2%in%(allpops[1:4]))] <- "resi_mig"
# Is between or within waterbody
Pops_long$Waterbody <- "inter"
Pops_long$Waterbody[substr(Pops_long$pop1, 1, 3)==substr(Pops_long$pop2, 1, 3)] <- "intra"

### Plot corrolation between Fst of all pop pairs
p <- ggplot(Pops_long[Pops_long$Waterbody=="intra",]) +
  geom_line(aes(start, Fst, group = Population, col = Population), alpha = 0.75) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
    geom_errorbar(data = Pops_long[1,],
    aes(xmin = 2e6, xmax = 22e6, y = -0.06), linewidth = 1, width = 0.05)+
  geom_text(data = Pops_long[1,],
    aes(x = ((22e6-2e6)/2)+2e6, y = -0.08), label = "20Mbp", size = 2) +
  theme_classic() +
  theme(legend.position = "none",panel.spacing = unit(0,'lines')) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  scale_color_manual(values = cbPalette) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

ggsave(paste0(my_bins, "_Fst_across_genome.png"),p,  width = 10, height=4)

### Plot corrolation between Fst of all pop pairs
p <- ggplot(Pops_long[Pops_long$Waterbody=="intra",]) +
  geom_line(aes(start, Fst, group = Population, col = pop1)) +
  facet_grid(pop1~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "none",panel.spacing = unit(0,'lines')) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  scale_color_manual(values = cbPalette, name = "Lagoon") +
  theme(panel.spacing = unit(0,'lines'), legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

ggsave(paste0(my_bins, "_Fst_across_genome_bypop.png"),p,  width = 10, height=5)
ggsave(paste0("test.png"),p,  width = 10, height=5)


