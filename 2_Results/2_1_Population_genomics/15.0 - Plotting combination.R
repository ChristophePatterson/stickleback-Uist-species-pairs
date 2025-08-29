# CalcuLation of LEA on stickleback popuLations
library(patchwork)
library(ggplot2)
library(ape)
library(vcfR)
library(tidyverse)
library(LEA)
library(ggnewscale)
# BiocManager::install(version = '3.20')
# BiocManager::install("LEA")

## Colorblind palette
cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

vcf.ver <- "GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20"
SNP.library.name <- "stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000"

# Test whether working on HPC or laptop and set working directory accordingly
# Laptop test
if(grepl(getwd(), pattern = "C:/Users/mbzcp2/")){
  dir.path <-"C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/"
  plot.dir <- "C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/"
}
# HPC test
if(grepl(getwd(), pattern = "/gpfs01/home/mbzcp2")){
  dir.path <-"/gpfs01/home/mbzcp2/data/sticklebacks/"
  plot.dir <- paste0("/gpfs01/home/mbzcp2/data/sticklebacks/results/", vcf.ver)
}
## Create directory is not already
dir.create(plot.dir)
dir.create(paste0(plot.dir, "/LEA_PCA/"))
dir.create(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/"))

## Read in vcf file
geno <- read.geno(paste0(dir.path, "results/",vcf.ver, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name, "_paired.geno"))
# Get an read sample information
samples_data <- data.frame(ID = readLines(paste0(dir.path, "results/",vcf.ver, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name, "_paired.geno.samples.txt")))
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)

dim(geno)

#########################
# # # #  MDS plot # # # #
#########################

## Run MDS
geno.mat <- geno
geno.mat[geno.mat==9] <- NA
dc <- dist(geno.mat)
mds <- cmdscale(dc, k = 4)

### Merge sample data with MDS
samples_data <- cbind.data.frame(samples_data, data.frame(mds))
colnames(samples_data) <- gsub("X", "MDS", colnames(samples_data))

## Rename anad to migratory
samples_data$Ecotype[samples_data$Ecotype=="anad"] <- "mig"

## MDS plots
mds12.plot <- ggplot(samples_data) +
  geom_point(aes(MDS1, MDS2, col = Waterbody, shape = Ecotype), size  = 3, show.legend = F) +
  scale_color_manual(values = cbPalette) +
  labs(x = "MDS1", y = "MDS2") + theme_bw() + theme(panel.border = element_rect(color = "black", linewidth = 1))

mds23.plot <- ggplot(samples_data) +
  geom_point(aes(MDS1, MDS3, col = Waterbody, shape = Ecotype), size  = 3) +
  scale_color_manual(values = cbPalette) +
  labs(x = "MDS2", y = "MDS3") + theme_bw() + theme(panel.border = element_rect(color = "black", linewidth = 1))
# Combine
mdsplot <- (mds12.plot + mds23.plot)

ggsave("test.png", mdsplot, width = 10, height = 5)

#########################
# # # #  PopHet  # # # #
#########################

my_bins <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/indPops/sliding_window_w100kb_s100kb_m1_PopPair_auto"

## Set order pops should be plotted
allpops <- c("CLAC" ,"DUIN","LUIB","OBSE",  "CLAM", "DUIM", "LUIM", "OBSM")

# Read in data
hetPops <- read_csv(paste0(my_bins,".csv"))

## Pivot pi
hetPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("pi_", colnames(hetPops))])], 
             cols = colnames(hetPops)[grep("pi_", colnames(hetPops))], values_to = "pi", names_prefix = "pi_", names_to = "Population")
## Pivot dxy
dxyPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("dxy_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("dxy_", colnames(hetPops))], values_to = "dxy", names_prefix = "dxy_", names_to = "Population")
## Pivot fst
FstPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("Fst_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("Fst_", colnames(hetPops))], values_to = "Fst", names_prefix = "Fst_", names_to = "Population")
## Pivot fst
HetInd_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("het_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("het_", colnames(hetPops))], values_to = "het", names_prefix = "het_", names_to = "individual")

## Summarise dxy
dxyPops_long_summary <- dxyPops_long %>%
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
  group_by(by = Population) %>%
  summarise(pi.mn = mean(pi, na.rm =T), pi.sd = sd(pi, na.rm = T), pi.md = median(pi, na.rm =T)) %>%
  rename(Pop = by)

hetPops_long_summary$Pop <- factor(hetPops_long_summary$Pop, levels = allpops)

## plot
p.pops <- ggplot() +
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
  theme_bw() +
  theme(panel.background = element_rect(fill = "grey"), legend.frame = element_rect(colour = 'black'),
        axis.title = element_blank(), axis.text.y = element_text(angle=90, hjust = 0.5)) +
    annotate("text", x = -0.25, y = 6.5, label = "Migratory", angle = 90) + ## Axis annotation
    annotate("text", x = -0.25, y = 2.5, label = "Resident", angle = 90) +
    annotate("text", y = -0.25, x = 6.5, label = "Migratory") +
    annotate("text", y = -0.25, x = 2.5, label = "Resident") +
    coord_fixed(clip = 'off', x = c(0.5, 8.5), y = c(0.5, 8.5))

#########################
# # # # #  Fst # # # # #
#########################


my_bins_fst <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi"
print(my_bins_fst)

# read in data
sliding_wd <- as_tibble(read.csv(paste0(my_bins_fst, ".csv"), header = T))

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

chr <- chr[chr$GenBank.seq.accession%in%unique(sliding_wd$scaffold),]
chr$Sequence.name <- gsub("chr", "", chr$Sequence.name)
# chr <- chr[order(chr$Seq.length, decreasing = T),]
chr$bi.col <- rep(1:2, length.out = dim(chr)[1])

sliding_wd$chr <- chr$Sequence.name[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]
sliding_wd$bi.col <- chr$bi.col[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]

## Rename Fst column to code is usualble across results
colnames(sliding_wd)[grep("Fst",colnames(sliding_wd))] <- "Fst"
colnames(sliding_wd)[grep("dxy",colnames(sliding_wd))] <- "dxy"

# Reorder facet labels
sliding_wd$chr <- factor(sliding_wd$chr, levels = chr$Sequence.name)

## Fst plot
p.fst <- ggplot(sliding_wd, aes(x = start, y = Fst, group = chr, col = as.factor(bi.col))) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = start, , col = as.factor(bi.col), fill = as.factor(bi.col)), linewidth = 0.2) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey60")) +
  scale_fill_manual(values = c("black", "grey60")) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  geom_errorbar(data = sliding_wd[sliding_wd$chr=="I",],
    aes(xmin = 2e6, xmax = 22e6, y = 0.7), linewidth = 1, width = 0.05)+
  geom_text(data = sliding_wd[sliding_wd$chr=="I",],
    aes(x = ((22e6-2e6)/2)+2e6, y = 0.75, label = "20Mbp"), size = 3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(sliding_wd$Fst, na.rm = T) + 0.05), sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

plot1 <- (mdsplot + p.pops)/p.fst + plot_layout(heights = c(10, 6)) + plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
## Save
ggsave(paste0("test.png"), plot1 , height = 9, width = 15.92)
ggsave(paste0(plot.dir, "/Figure_1.png"),plot1 , height = 9, width = 15.92)
ggsave(paste0(plot.dir, "/Figure_1.pdf"),plot1 , height = 9, width = 15.92)
