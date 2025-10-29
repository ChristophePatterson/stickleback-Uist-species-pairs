## Read in bed files and plot which regions were not lifted over

## load libraries
library(tidyverse)
library(patchwork)

## Read in scaffold info to get lengths for plotting

## Chromosome levels
chr_levels <- c(as.character(as.roman(c(1:21))),"Un")

## For DUKE genome
chr_DUKE <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
## Calculate cumulative position
chr_DUKE$Cum.Seq.length <- c(0, cumsum(chr_DUKE$Seq.length[1:(nrow(chr_DUKE)-1)]))

## For gasAcu5 genome
chr_v5 <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T))

# Sort out factor order of chromosomes
chr_v5 <- chr_v5 %>%
    filter(Chromosome.name %in% c(as.character(as.roman(c(1:21))), "Un")) %>%
    mutate(chr = factor(Chromosome.name, levels=chr_levels))

## Calculate cumulative position
chr_v5$Cum.Seq.length <- c(0, cumsum(chr_v5$Seq.length[1:(nrow(chr_v5)-1)]))

## Adjust length of "Un" chromosome to end at total genome length
chr_v5 <- chr_v5 %>%
    mutate(Seq.length = ifelse(Chromosome.name=="Un", Cum.Seq.length - sum(Seq.length[Chromosome.name!="Un"]), Seq.length))

# Sort out factor order of chromosomes
chr_DUKE <- chr_DUKE %>%
    filter(Chromosome.name %in% c(1:21)) %>%
    mutate(chr = factor(as.character(as.roman(Chromosome.name)), levels=chr_levels))

## Read in bed files for original regions and lifted over regions

## Roberts et al EcoPeaks
EcoPeaks_v5 <- read.table("Roberts-et-al-2021-Specific-EcoPeaks-v5.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, ID=V4) %>%
    mutate(chr = factor(gsub("chr", "", chr), levels=chr_levels),
           data = "Roberts et al EcoPeaks")
## Read in regions that were lifted over to DUKE
EcoPeaks_DUKE <- read.table("Roberts-et-al-2021-Specific-EcoPeaks-vDUKE.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, ID=V4) %>%
    mutate(chr = factor(as.character(as.roman(gsub("chr|chr0", "", chr))), levels=chr_levels),
           data = "Roberts et al EcoPeaks")

## Jones et al CSS
Jones_v5 <- read.table("Jones-et-al-2012-CSS-02-v4_wID.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, ID=V4) %>%
    mutate(chr = factor(gsub("chr|chr0", "", chr), levels=chr_levels),
    data = "Jones et al CSS")
## Read in regions that were lifted over to DUKE
Jones_DUKE <- read.table("Jones-et-al-2012-CSS-02-vDUKE.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, ID=V4) %>%
    mutate(chr = factor(as.character(as.roman(gsub("chr|chr0", "", chr))), levels=chr_levels),
           data = "Jones et al CSS")

## Read in v5 gene annotations for plotting translocation positions
genes_v5 <- read.table("../stickleback_v5_ensembl_genes.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, type=V4, ID=V5) %>%
    mutate(chr = factor(gsub("chr", "", chr), levels=chr_levels),
           ID = gsub("ID=", "", str_split_i(ID, ";", 1)),
           data = "Genes v5")

## Read in DUKE gene annotations for plotting translocation positions
genes_DUKE <- read.table("../stickleback_DUKE_ensembl_genes.bed", header=F) %>%
    rename(chr = V1, start = V2, end = V3, type=V4, ID=V7) %>%
    mutate(chr = factor(as.character(as.roman(gsub("chr|chr0", "", chr))), levels=chr_levels),
           ID = gsub("ID=", "", str_split_i(ID, ";", 1)),
           data = "Genes DUKE")

## Which regions were not lifted over
EcoPeaks_v5$unlifted <- ifelse(EcoPeaks_v5$ID%in%EcoPeaks_DUKE$ID, "Lifted", "Not_lifted")
Jones_v5$unlifted <- ifelse(Jones_v5$ID %in% Jones_DUKE$ID, "Lifted", "Not_lifted")
genes_v5$unlifted <- ifelse(genes_v5$ID %in% genes_DUKE$ID, "Lifted", "Not_lifted")

## Combine data for plotting
p <- ggplot() +
    geom_segment(data=chr_v5, aes(y=Chromosome.name, yend=Chromosome.name, x=0, xend=Seq.length), color="grey60", linewidth=5, lineend = "round") +
    geom_segment(data = EcoPeaks_v5, aes(y=chr, yend = chr, x=start, xend = end, color=unlifted), lineend = "round", linewidth = 5, alpha = 0.5) +
    geom_segment(data = Jones_v5, aes(y=chr, yend = chr, x=start, xend = end, color=unlifted), lineend = "round", linewidth = 5, alpha = 0.5) +
    theme_bw() +
    facet_wrap(~data, ncol=1) +
    labs(
         y="Chromosome",
         x="Start position (bp)") +
    scale_color_manual(values=c("Lifted"="grey10", "Not_lifted"="red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("EcoPeaks_Jones_liftOver_success_summary.png", plot=p, width=12, height=8)

## Combine data for plotting
p <- ggplot() +
    geom_segment(data=chr_v5, aes(y=Chromosome.name, yend=Chromosome.name, x=0, xend=Seq.length), color="grey60", linewidth=5, lineend = "round") +
    geom_segment(data = genes_v5[genes_v5$type=="gene",], aes(y=chr, yend = chr, x=start, xend = end, color=unlifted), linewidth = 5, alpha = 0.5) +
    theme_bw() +
    labs(
         y="Chromosome",
         x="Start position (bp)") +
    scale_color_manual(values=c("Lifted"="grey10", "Not_lifted"="red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("EcoPeaks_Jones_liftOver_gene_success_summary.png", plot=p, width=12, height=8)

## Merge locations of lifted over regions from v4 and v5
EcoPeaks <- merge(EcoPeaks_v5, EcoPeaks_DUKE, by="ID") %>%
    rename(chr_v5 = chr.x, start_v5 = start.x, end_v5 = end.x,
           chr_DUKE = chr.y, start_DUKE = start.y, end_DUKE = end.y,
           chr = chr.x)

# Plot lifted over regions
p.EcoPeaks <- ggplot(EcoPeaks) +
    geom_segment(data=chr_v5, aes(y="v5", yend="v5", x=0, xend=Seq.length), color="grey60", linewidth=2, lineend = "round") +
    geom_segment(data=chr_DUKE, aes(y="DUKE", yend="DUKE", x=0, xend=Seq.length), color="grey60", linewidth=2, lineend = "round") +
    geom_segment(aes(y = "v5", yend = "DUKE", x = start_v5, xend = start_DUKE), alpha = .5)+
    geom_segment(aes(y = "v5", yend = "DUKE", x = end_v5, xend = end_DUKE), alpha = .5)+
    geom_segment(data = EcoPeaks_v5, aes(y = "v5", yend = "v5",, x=start, xend = end, color=unlifted), lineend = "round", linewidth = 2) +
    geom_segment(data = EcoPeaks_DUKE, aes(y = "DUKE", yend = "DUKE",, x=start, xend = end), color="grey10", lineend = "round", linewidth = 2) +
    theme_bw() +
    facet_grid(chr~., drop = FALSE, , scales = "free_x") +
    scale_color_manual(values=c("Lifted"="grey10", "Not_lifted"="red")) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6), name = "Mbs") +
    labs(col = "Success of LiftOver",
         y="Chromosome",
         x="Start position (bp)") +
    theme(panel.spacing = unit(0,'lines'), legend.position = "bottom",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1)) +
    ggtitle("Roberts et al 2021: Specific EcoPeaks")

## Merge locations of lifted over regions from v4 and v5
Jones <- merge(Jones_v5, Jones_DUKE, by="ID") %>%
    rename(chr_v5 = chr.x, start_v5 = start.x, end_v5 = end.x,
           chr_DUKE = chr.y, start_DUKE = start.y, end_DUKE = end.y,
           chr = chr.x)

# Plot lifted over regions
p.jones <- ggplot(Jones) +
    geom_segment(data=chr_v5, aes(y="v5", yend="v5", x=0, xend=Seq.length), color="grey60", linewidth=2, lineend = "round") +
    geom_segment(data=chr_DUKE, aes(y="DUKE", yend="DUKE", x=0, xend=Seq.length), color="grey60", linewidth=2, lineend = "round") +
    geom_segment(aes(y = "v5", yend = "DUKE", x = start_v5, xend = start_DUKE), alpha = .5)+
    geom_segment(aes(y = "v5", yend = "DUKE", x = end_v5, xend = end_DUKE), alpha = .5)+
    geom_segment(data = Jones_v5, aes(y = "v5", yend = "v5",, x=start, xend = end, color=unlifted), lineend = "round", linewidth = 2) +
    geom_segment(data = Jones_DUKE, aes(y = "DUKE", yend = "DUKE",, x=start, xend = end), color="grey10", lineend = "round", linewidth = 2) +
    theme_bw() +
    facet_grid(chr~., drop = FALSE, scales = "free_x") +
    scale_color_manual(values=c("Lifted"="grey10", "Not_lifted"="red")) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6), name = "Mbs") +
    labs(col = "Success of LiftOver",
         y="Chromosome",
         x="Start position (bp)") +
    theme(panel.spacing = unit(0,'lines'), legend.position = "bottom",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1)) +
    ggtitle("Jones et al 2012: CSS_02")


ggsave("EcoPeaks_Jones_liftOver_success_with_translocation_poistion.png",
    plot= ((p.EcoPeaks + p.jones) / guide_area()) + plot_layout(heights = c(30,1), guides = "collect"), 
    height = 24.62*0.75, width = 7.96*2)
