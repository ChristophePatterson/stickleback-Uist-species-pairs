# CalcuLation of LEA on stickleback popuLations
library(patchwork)
library(ggplot2)
library(ape)
library(vcfR)
library(tidyverse)
library(LEA)
library(ggnewscale)
library(ggrepel)
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
        axis.title = element_blank(), axis.text.y = element_text(angle=90, hjust = 0.5),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
    annotate("text", x = -0.25, y = 6.5, label = "Migratory", angle = 90) + ## Axis annotation
    annotate("text", x = -0.25, y = 2.5, label = "Resident", angle = 90) +
    annotate("text", y = -0.25, x = 6.5, label = "Migratory") +
    annotate("text", y = -0.25, x = 2.5, label = "Resident") +
    coord_fixed(clip = 'off', x = c(0.5, 8.5), y = c(0.5, 8.5)) +
    guides(color = guide_legend(override.aes = list(size = 0.2)))

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

chr$Sequence.name <- factor(chr$Sequence.name, levels = chr$Sequence.name)
sliding_wd$chr <- chr$Sequence.name[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]
sliding_wd$bi.col <- chr$bi.col[match(sliding_wd$scaffold, chr$GenBank.seq.accession)]

chr <- chr[order(chr$Sequence.name),]
# Calculate cumulative length (starting at 0 )
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))
# Amened onto CSS data
sliding_wd$start.cum <- sliding_wd$start+(chr$Cum.Seq.length[match(sliding_wd$chr, factor(chr$Sequence.name, levels = levels(sliding_wd$chr)))]+chr$Cum.Seq.length[1])

## Rename Fst column to code is usualble across results
colnames(sliding_wd)[grep("Fst",colnames(sliding_wd))] <- "Fst"
colnames(sliding_wd)[grep("dxy",colnames(sliding_wd))] <- "dxy"

# Reorder facet labels
sliding_wd$chr <- factor(sliding_wd$chr, levels = chr$Sequence.name)

## Fst plot
p.fst <- ggplot(sliding_wd) +
  geom_vline(xintercept = chr$Cum.Seq.length, col = "grey80") +
  geom_text(data = chr, aes(x = Cum.Seq.length+(Seq.length/2), y = 0.95, label = Sequence.name)) +
  #geom_point(show.legend = F, size = 0.2) +
  geom_ribbon(aes(ymax = Fst, ymin = 0, x = start.cum, , col = as.factor(bi.col), fill = as.factor(bi.col), group = chr), linewidth = 0.2, show.legend = F) +
  # geom_line(aes(x = mid, y = rolling_av-1.05), col = "red", size = 1.5, linewidth = 4) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_fill_manual(values = c("black", "grey50")) +
  # facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  # geom_errorbar(data = sliding_wd[sliding_wd$chr=="I",],
  #   aes(xmin = 2e6, xmax = 22e6, y = 0.7), linewidth = 1, width = 0.05)+
  # geom_text(data = sliding_wd[sliding_wd$chr=="I",],
  #   aes(x = ((22e6-2e6)/2)+2e6, y = 0.75, label = "20Mbp"), size = 3) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(sliding_wd$start.cum),20e6)),name = "Mbps", expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  ylab("Fst") +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1))

        
plot1 <- (mdsplot + p.pops)/p.fst + plot_layout(heights = c(10, 6)) + plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")
## Save
ggsave(paste0("test.png"), plot1 , height = 8, width = 15.92)
ggsave(paste0(plot.dir, "/Figure_1.png"),plot1 , height = 9, width = 15.92)
ggsave(paste0(plot.dir, "/Figure_1.pdf"),plot1 , height = 9, width = 15.92)


##############################
 # # # # # # CSS # # # # # # 
##############################

## Read in combined CSS scores
CSS.dir <-  "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/stickleback.wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05"
CSS.run <- "stickleback.wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05.CSSm.10000perm.txt"

CSS <- read.table(paste0(CSS.dir,"/", CSS.run),header = T)

# Add chr info to CSS
CSS$chr <- chr$Sequence.name[match(CSS$chr, chr$GenBank.seq.accession)]
CSS$bi.col <- chr$bi.col[match(CSS$chr, chr$Sequence.name)]

# Add cumulative sum to CSS
CSS$start.cum <- CSS$start+(chr$Cum.Seq.length[match(CSS$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])

# Filter out windows with low or high number of SNPs
CSS.HQ <- CSS[CSS$nsnps>=5&CSS$nsnps<=200,]

##  Read in dropped populations
CSS.drop <- read_csv("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_combine.csv")
# Get regions that are signicant across all regions
CSS.drop.all.sig <- na.omit(CSS.drop[CSS.drop$all.sig.qvalue.0001,])

## Combine CSS and CSS dropPop
CSS.HQ$wnd.name <- paste(CSS.HQ$chr, CSS.HQ$start, CSS.HQ$end, sep =  "_")
CSS.drop.all.sig$wnd.name <- paste(CSS.drop.all.sig$chr, CSS.drop.all.sig$start, CSS.drop.all.sig$end, sep =  "_")

## Which regions are held sig across all droppops
CSS.HQ$drop.all.sig.qvalue.0001 <- CSS.HQ$wnd.name %in% CSS.drop.all.sig$wnd.name

## Read in CSS regions with gene annotations
CSS.annotations <- read_csv("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions_grouped.txt")
CSS.annotations$chr <- factor(CSS.annotations$chr, levels = levels(CSS.HQ$chr))

## Which are the highest CSS sig regions
CSS.annotations.top.regions <- CSS.annotations# [CSS.annotations$mn.CSS>=2,]

##  Select region of interest
# regions <- data.frame(chr = factor(c("I", "IV", "XI", "XXI"), levels = levels(chr$Sequence.name)), start = c(25000000, 12000000, 5000000, 8000000), end = c(31000000, 16000000, 10000000, 15000000))
regions <- data.frame(chr = factor(c("I", "IV", "XIX", "XXI"), levels = levels(chr$Sequence.name)), start = c(26000000, 12000000, 2000000, 8000000), end = c(28000000, 16000000, 10000000, 15000000))
regions$start.cum <- regions$start+(chr$Cum.Seq.length[match(regions$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])
regions$end.cum <- regions$end+(chr$Cum.Seq.length[match(regions$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])

# Filer dataset to specific region
CSS.HQ.filt <- CSS.HQ %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

# Filer dataset to specific region
CSS.annotations.top.regions.filt <- CSS.annotations.top.regions %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

# Plot
p.CSS <- ggplot(CSS.HQ[!CSS.HQ$drop.all.sig.qvalue.0001,]) +
  geom_vline(xintercept = chr$Cum.Seq.length, col = "grey80") +
  geom_segment(data = regions, aes(x = start.cum, xend = end.cum, y = -max(CSS.HQ$css)*0.05), linewidth = 3) +
  geom_text(data = chr, aes(x = Cum.Seq.length+(Seq.length/2), y = max(CSS.HQ$css)*1.1, label = Sequence.name)) +
  geom_point(aes(start.cum, css, col = as.factor(bi.col)), show.legend = F) +
  geom_point(data = CSS.HQ[CSS.HQ$drop.all.sig.qvalue.0001,], aes(start.cum, css), col = "firebrick3") +
  scale_color_manual(values = c("black", "grey50")) +
  scale_fill_manual(values = c("black", "grey50")) +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(sliding_wd$start.cum),20e6)),name = "Mbps", expand = c(0,0)) +
  ylab("CSS") +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1))

# Save
CSS.annotations.top.regions.filt$genes <- gsub("\\|", ",", CSS.annotations.top.regions.filt$genes)
CSS.annotations.top.regions.filt$genes.filt <- gsub("NA", "",do.call("c", lapply(str_split(CSS.annotations.top.regions.filt$genes, ","),
                                                 function(x) paste(x[!grepl("ENSGACG", x)], collapse = ", "))))


p.CSS.filt <- ggplot(CSS.HQ.filt[!CSS.HQ.filt $drop.all.sig.qvalue.0001,]) +
  geom_point(aes(start, css, col = as.factor(bi.col)), show.legend = F) +
  geom_point(data = CSS.HQ.filt [CSS.HQ.filt $drop.all.sig.qvalue.0001,], aes(start, css), col = "firebrick3") +
  geom_text_repel(data = CSS.annotations.top.regions.filt, 
          aes(x = start+((end-start)/2), y = mn.CSS, label = gsub(", ", "\n", genes.filt)), hjust = 0, nudge_y = 5, nudge_x = 250000,
          direction = "both", box.padding = 0.1,
          size = 1.5, max.overlaps = 15, min.segment.length = 0) +
  geom_segment(data = CSS.annotations.top.regions.filt, 
        aes(x = start, xend = end, y = -1)) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_fill_manual(values = c("black", "grey50")) +
  facet_wrap(~chr, scale = "free_x") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(sliding_wd$start),0.5e6)),name = "Mbps", expand = c(0,0)) +
  ylab("CSS") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1))

CSS.plot.comb <- p.CSS/p.CSS.filt + plot_layout(heights=c(1,2)) + plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")")


# ggsave("test.png", CSS.plot.comb , height = 15.92*0.66666, width = 15.92*0.66666)
ggsave(paste0(plot.dir, "/Figure_CSS.pdf"), CSS.plot.comb , height = 15.92*0.66666, width = 15.92*0.66666)
ggsave(paste0(plot.dir, "/Figure_CSS.png"), CSS.plot.comb , height = 15.92*0.66666, width = 15.92*0.66666)

