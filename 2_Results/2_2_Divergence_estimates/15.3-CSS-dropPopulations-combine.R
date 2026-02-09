## devtools::install_github("jdstorey/qvalue")

## Run this code as
# module load R-uoneasy/4.2.1-foss-2022a
# Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/15.3-CSS-dropPopulations-combine.R &> /gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops/DropPops_combined_log.txt

library(qvalue)
library(tidyverse)
library(patchwork)
library(IRanges)

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

## args<-commandArgs(trailingOnly=T)
## CSS.dir <- args[1]
CSS.dir <-  "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops"

CSS.run <- ".wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05"

## Get combinations of dropped pops
dropComb <- data.frame(waterbody = c("CLAC", "DUIN", "LUIB", "OBSE"))

dropComb$dropped <- paste0("n", dropComb$waterbody)
# Getfile locations
dropComb$run_name <- paste0("stickleback.", dropComb$dropped, CSS.run)
dropComb$run_dir <- paste0(CSS.dir, "/", dropComb$dropped, "/", dropComb$run_name, "/")
dropComb$run_file <- paste0(dropComb$run_dir, dropComb$run_name, ".CSSm.10000perm.txt")
dropComb$run_file

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
chr$Sequence.name <- gsub("chr", "", chr$Sequence.name)

# Calculate cumumative length
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))
## Read in each CSS calc
CSS.list <- apply(dropComb, 1, function(x) {
    df <- cbind.data.frame(
      dropped = x["dropped"],
      read.table(x["run_file"], header = TRUE)
    )
    df <- df[df$nsnps>=5|df$nsnps<=200,]
    df$qval.0001 <- qvalue(1 - df$pval, fdr.level = 0.0001)$qvalues
    df
  })

# Combined into single data frame
CSS.long <- do.call("rbind", CSS.list)

CSS.long$chr <- chr$Sequence.name[match(CSS.long$chr, chr$GenBank.seq.accession)]
CSS.long$chr <- factor(CSS.long$chr, levels = chr$Sequence.name)

CSS.long <- CSS.long %>%
    mutate(qval.sig.0001 = qval.0001 < 0.001)

CSS.long <- CSS.long %>%
  group_by(dropped) %>%
  mutate(temp = replace(qval.sig.0001, is.na(qval.sig.0001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1, .add = T) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  group_by(dropped, .add = T) %>%
  mutate(goal.0001 = ifelse(temp, paste0(dplyr::first(dropped), "_", dplyr::first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

top.regions <- sort(table(CSS.long$goal.0001[CSS.long$css>=2]), decreasing = T)*500

## Write out combined CSS calculations
head(CSS.long)
CSS.long  %>% 
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_combine.csv"), row.names = F, quote = F, sep = ",")

print("Table of each qval for each dropped population")
table(CSS.long$dropped, CSS.long$qval.sig.0001, useNA="ifany")

# Pivot table wider
CSS.wide <- CSS.long %>%
  pivot_wider(
    names_from = dropped,
    values_from = c(nsnps, css, nperms, pval, qval.sig.0001, qval.0001, goal.0001)
  )

## png(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_cor.png"), width = 1000, height = 1000)
## plot(CSS.wide[,grepl("css", colnames(CSS.wide))])
## dev.off()

p <- ggplot(CSS.long) +
  geom_point(aes(start, css, col = qval.sig.0001), size = 0.5) +
  scale_color_manual(values=c("grey50", "deepskyblue"), name = "qvalue\nsignificance\n(<0.0001)") +
  geom_point(data = CSS.long[CSS.long$goal.0001%in%names(top.regions[top.regions>=2500]),], aes(start, css), col = "red", size = 0.5) +
  facet_grid(dropped~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(10e6, max(chr$Seq.length),10e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "bottom",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

# ggsave("test.png", p, width = 7.96, height = 7.96/1.8)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS.png"), p, , width = 7.96*2, height = 7.96)

## p <- ggplot(CSS.long) +
##     geom_segment(aes(x = start, xend = end, y = dropped, yend= dropped, col = qval.sig.0001), linewidth = 3) +
##     geom_segment(data = CSS.long[CSS.long$goal.0001%in%names(top.regions[top.regions>=2500]),], aes(x = start, xend = end, y = dropped, yend= dropped), col = "red", linewidth = 3) +
##     facet_grid(chr~., scale = "free_x")
## 
## ggsave("test.png", p, width = 8, height = 15)
regions <- data.frame(chr = c("I", "IV", "XI", "XXI"), start = c(25000000, 12000000, 5000000, 8000000), end = c(31000000, 16000000, 10000000, 15000000))

regions <- data.frame(chr = factor(c("I", "IV", "IX", "XI", "XXI"), levels = chr$Sequence.name), start = c(26000000, 12000000, 12000000, 5500000, 9000000), end = c(27500000, 16000000, 14500000,7000000, 12500000))

# Filer dataset to specific region
CSS.long.filt <- CSS.long %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

regions_plot <- ggplot(CSS.long.filt) +
  geom_point(aes(start, css, col = qval.sig.0001), size = 0.5) +
  scale_color_manual(values=c("grey50", "firebrick3"), name = "qvalue significance\n(<0.0001)") +
  facet_grid(dropped~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),1e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank(), legend.position = "bottom")

# ggsave("test.png", regions_plot, width = 7.96*2, height = 7.96)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_regions.png"), regions_plot, width = 7.96*2, height = 7.96)


#### CALULATE REGIONS THAT ARE SIGNIFICANT ACCROSS ALL dropPops
# Pivot table wider
CSS.wide <- CSS.long %>%
  pivot_wider(
    names_from = dropped,
    values_from = c(goal, nsnps, css, nperms, pval, qval.sig.0001, qval.0001, goal.0001)
  )

## png(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_cor.png"), width = 1000, height = 1000)
## plot(CSS.wide[,grepl("css", colnames(CSS.wide))])
## dev.off()

## Determine which windows are significant across all dropped population comparisons
CSS.wide$all.sig.qvalue.0001 <- apply(CSS.wide[,grepl("qval.sig.0001_", colnames(CSS.wide))], MARGIN = 1, function(x) 
  ifelse(any(is.na(x)), return(NA), all(x)))

## Determine which windows are significant across any dropped population comparisons
CSS.wide$any.sig.qvalue.0001 <- apply(CSS.wide[,grepl("qval.sig.0001_", colnames(CSS.wide))], MARGIN = 1, function(x) any(x))

print("XXXXXXXXX STATS XXXXXXXXXXXX")
print("Data proportions")
dim(CSS.wide[,grepl("qval.sig.0001_", colnames(CSS.wide))])
print("Total number of sig windows across all genome for all dropped pops")
table(CSS.wide$all.sig.qvalue.0001, useNA="ifany")
print("Total number of sig windows that were sig for at least one drop pop")
table(CSS.wide$any.sig.qvalue.0001, useNA="ifany")

# Pivot a subset of this to be in tidy format
CSS.long.all.sig <- pivot_longer(CSS.wide[,c("chr","start","end", paste0("css_", dropComb$dropped), "all.sig.qvalue.0001")],
            cols = paste0("css_", dropComb$dropped), values_to = "css", names_prefix = "css_", names_to = "dropped")
# Add on individual qvalues
CSS.long.all.sig <- cbind(CSS.long.all.sig, pivot_longer(CSS.wide[,c(paste0("qval.0001_", dropComb$dropped))],
            cols = paste0("qval.0001_", dropComb$dropped), values_to = "qval.0001", names_prefix = "qval.0001_", names_to = "dropped")[,2])

CSS.long.all.sig  %>% 
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_combine.csv"), row.names = F, quote = F, sep = ",")

# Plot
p.all <- ggplot(na.omit(CSS.long.all.sig[!CSS.long.all.sig$all.sig.qvalue.0001,])) +
  geom_point(aes(start, css, col = qval.0001<0.0001), size = 0.5) +
  geom_point(data = na.omit(CSS.long.all.sig[CSS.long.all.sig$all.sig.qvalue.0001,]), aes(start, css), size = 0.5, col = "red") +
  scale_color_manual(values=c("grey50", "deepskyblue"), name = "qvalue\nsignificance\n(<0.0001)") +
  facet_grid(dropped~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(10e6, max(chr$Seq.length),10e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "bottom",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig.png"), p, width = 24.62/2, height = 7.96)
# ggsave("test.png", p.all, width = 24.62/2, height = 7.96)

# Filer dataset to specific region
CSS.long.all.sig.filt <- CSS.long.all.sig %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

# Plot
p.all.regions <- ggplot(na.omit(CSS.long.all.sig.filt[!CSS.long.all.sig.filt$all.sig.qvalue.0001,])) +
  geom_point(aes(start, css, col = qval.0001<0.0001), size = 0.5) +
  geom_point(data = na.omit(CSS.long.all.sig.filt[CSS.long.all.sig.filt$all.sig.qvalue.0001,]), aes(start, css), size = 0.5, col = "red") +
  scale_color_manual(values=c("grey50", "deepskyblue"), name = "qvalue\nsignificance\n(<0.0001)") +
  facet_grid(dropped~chr, scale = "free_x", space = "free_x") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),2e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

# ggsave("test.png", p.all.regions, width = 7.96*2, height = 7.96)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig.regions.png"), p.all.regions, width = 7.96*2, height = 7.96)

## Find all the contigous regions
CSS.wide <- CSS.wide %>%
  mutate(temp = replace(all.sig.qvalue.0001, is.na(all.sig.qvalue.0001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1, .add = T) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  mutate(all.goal.0001 = ifelse(temp, paste0(dplyr::first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)


# Calculate mean across all populations
CSS.wide$mn.CSS <- apply(CSS.wide, MARGIN = 1, function(x) mean(as.numeric(x[grep("css_", names(x))]), na.rm = T))
CSS.wide$mx.CSS <- apply(CSS.wide, MARGIN = 1, function(x) max(as.numeric(x[grep("css_", names(x))]), na.rm = T))

## Table
top.regions.table <- CSS.wide %>%
  filter(!is.na(all.goal.0001)) %>%
  group_by(all.goal.0001) %>%
  summarise(chr = dplyr::first(chr),
            start = min(start),
            end = max(end),
            mn.CSS = mean(mn.CSS),
            mx.CSS = max(mx.CSS),
            wnd.length = (1+max(end))-min(start)) %>%
  ungroup() %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(bp.break = c(Inf, diff(start)))

#### Print stats
print("############ STATS ##############")
print("Total contiguous regions")
length(unique(top.regions.table$all.goal.0001))
print("Total contiguous regions that are longer than 2500bp")
length(unique(top.regions.table$all.goal.0001[top.regions.table$wnd.length>2500]))
print("Total contiguous regions that are longer than 4500bp")
length(unique(top.regions.table$all.goal.0001[top.regions.table$wnd.length>4500]))
print("Total contiguous regions with mean CSS greater than 2")
length(unique(top.regions.table$all.goal.0001[top.regions.table$mn.CSS>2]))
print("Total contiguous regions with a CSS that exceeds 2")
length(unique(top.regions.table$all.goal.0001[top.regions.table$mx.CSS>2]))

# Write out table
top.regions.table %>% 
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_top_regions.txt"), row.names = F, quote = F, sep = ",")

# Remove regions that are smaller than 2500bp (less than three contigous significant windows)
top.regions.all <- sort(table(CSS.wide$all.goal.0001), decreasing = T)*500
top.regions.all <- top.regions.all[top.regions.all>=2500]

CSS.wide$top.regions <- CSS.wide$all.goal.0001 %in% names(top.regions.all)

## Find all the genes that are with these regions of high coverage
# Load in gene data
DUKE.bed <- tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/stickleback_DUKE_ensembl_genes.bed", header = F)) %>%
  dplyr::rename(chr = V1, start = V2, end = V3,type = V4, splitID=V5, splitID2=V6, ID = V7) %>%
  mutate(chr = factor(as.character(as.roman(as.numeric(gsub("chr", "", chr)))), levels = levels(CSS.long$chr)),
         gene.name = gsub("Name=|Parent=", "", str_split_i(ID, ";", 2)))

# Subset to genes
DUKE.bed.genes <- DUKE.bed[DUKE.bed$type == "gene",]

## Create  ranges overlap for all populations
DUKE.bed.ranges <- split(IRanges(DUKE.bed.genes$start, DUKE.bed.genes$end), DUKE.bed.genes$chr)

hits.tmp <- queryHits(findOverlaps(split(IRanges(top.regions.table$start, top.regions.table$end), top.regions.table$chr),
split(IRanges(DUKE.bed.genes[2271,]$start, DUKE.bed.genes[2271,]$end), DUKE.bed.genes[2271,]$chr)))

# Create blank columns
top.regions.table$contains.genes <- ""
top.regions.table$contains.genes.in10Kbps <- ""
top.regions.table$contains.genes.in100Kbps <- ""

# Loop through each region and ask which genes overlap with it
i <- 51
for(i in 1:nrow(top.regions.table)){
  # Ranges object for region
  ## Clear results from previous loop
  genes.tmp <- ""
  contains.genes.in10Kbps <- ""
  contains.genes.in100Kbps <- ""

  # Extract chromosome
  chr.tmp <- top.regions.table$chr[i]
  top.regions.range.tmp <- split(IRanges(top.regions.table$start[i], top.regions.table$end[i]), chr.tmp)
  # Overlaps
  DUKE.bed.hits.tmp <- findOverlaps(top.regions.range.tmp, DUKE.bed.ranges)[[chr.tmp]]
  # Extract gene names of overlaps
  genes.tmp <- ((DUKE.bed.genes[DUKE.bed.genes$chr==chr.tmp,])[subjectHits(DUKE.bed.hits.tmp),])$gene.name
  # Remove duplicates
  genes.tmp <- unique(genes.tmp)
  
  # Create flanking regions
  top.regions.flanks.tmp <- merge(flank(top.regions.range.tmp, start = TRUE, width = 10000), flank(top.regions.range.tmp, start = FALSE, width = 10000)) 
  DUKE.bed.hits.tmp <- findOverlaps(DUKE.bed.ranges, top.regions.flanks.tmp)[[chr.tmp]]
  genes.tmp.10Kbps <- (DUKE.bed.genes[queryHits(DUKE.bed.hits.tmp),])$gene.name
  # Remove genes already idenfied
  genes.tmp.10Kbps <- genes.tmp.10Kbps[!genes.tmp.10Kbps%in%genes.tmp]
  
  # Create flanking regions for 100000
  top.regions.flanks.tmp <- merge(flank(top.regions.range.tmp, start = TRUE, width = 100000), flank(top.regions.range.tmp, start = FALSE, width = 100000)) 
  DUKE.bed.hits.tmp <- findOverlaps(DUKE.bed.ranges, top.regions.flanks.tmp)[[chr.tmp]]
  genes.tmp.100Kbps <- (DUKE.bed.genes[queryHits(DUKE.bed.hits.tmp),])$gene.name
  # Remove genes already idenfied
  genes.tmp.100Kbps <- genes.tmp.100Kbps[!(genes.tmp.100Kbps%in%genes.tmp|genes.tmp.100Kbps%in%genes.tmp.10Kbps)]
  # Add to dataframe
  top.regions.table$contains.genes[i] <- paste0(genes.tmp, collapse = "|")
  top.regions.table$contains.genes.in10Kbps[i] <- paste0(genes.tmp.10Kbps, collapse = "|")
  top.regions.table$contains.genes.in100Kbps[i] <- paste0(genes.tmp.100Kbps, collapse = "|")
  # print(i)
}

# Write out file
top.regions.table %>% 
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_top_regions.txt"), row.names = F, quote = F, sep = ",")

#### Merge into segement separated by less that 100Kbps
# Create ranges object
top.regions.table.ranges <- split(IRanges(start = top.regions.table$start, end = top.regions.table$end), top.regions.table$chr)
    
# Merge with maxgap = 100000
top.regions.table.redueced <- reduce(top.regions.table.ranges, min.gapwidth = 100000)

# Find rows that were merged
top.regions.table.hits <- findOverlaps(top.regions.table.ranges, top.regions.table.redueced)

top.regions.table.hits.start <- start(top.regions.table.redueced)
top.regions.table.hits.end <- end(top.regions.table.redueced)

## Merge into data set  
top.grouped.regions.table <- data.frame(start = do.call("c", start(top.regions.table.redueced)), end = do.call("c", end(top.regions.table.redueced)))
top.grouped.regions.table$group <- rownames(top.grouped.regions.table)
top.grouped.regions.table$chr <- gsub('[0-9]', '', top.grouped.regions.table$group)

top.grouped.regions.table.genehits <- findOverlaps(top.regions.table.redueced, DUKE.bed.ranges)

top.grouped.regions.table$genes <- ""

for(i in 1:nrow(top.grouped.regions.table)){
  # Ranges object for region
  ## Clear results from previous loop
  genes.tmp <- ""

  # Extract chromosome
  chr.tmp <- top.grouped.regions.table$chr[i]
  top.regions.range.tmp <- split(IRanges(top.grouped.regions.table$start[i], top.grouped.regions.table$end[i]), chr.tmp)
  # Overlaps
  DUKE.bed.hits.tmp <- findOverlaps(top.regions.range.tmp, DUKE.bed.ranges)[[chr.tmp]]
  # Extract gene names of overlaps
  genes.tmp <- ((DUKE.bed.genes[DUKE.bed.genes$chr==chr.tmp,])[subjectHits(DUKE.bed.hits.tmp),])$gene.name
  # remove duplicates
  genes.tmp <- unique(genes.tmp)
  top.grouped.regions.table$genes[i] <- paste(genes.tmp, collapse = "|")

}

# Rearragne data set
top.grouped.regions.table$mn.CSS <- NA
top.grouped.regions.table$mn.CSS.sig <- NA
# add in CSS info
for(i in 1:nrow(top.grouped.regions.table)){
  start.tmp <- top.grouped.regions.table$start[i]
  end.tmp <- top.grouped.regions.table$end[i]
  chr.tmp <- top.grouped.regions.table$chr[i]
  CSS.tmp <- CSS.long[CSS.long$chr==chr.tmp&(between(CSS.long$start, start.tmp, end.tmp)|between(CSS.long$end, start.tmp, end.tmp)),]
  top.grouped.regions.table$mn.CSS[i] <- mean(CSS.tmp$css)
  top.grouped.regions.table$mn.CSS.sig[i] <- mean(CSS.tmp$css[CSS.tmp$qval.sig.0001])
}

# Add cumulative position
top.grouped.regions.table$start.cum <- top.grouped.regions.table$start+(chr$Cum.Seq.length[match(top.grouped.regions.table$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])
top.grouped.regions.table$end.cum <- top.grouped.regions.table$end+(chr$Cum.Seq.length[match(top.grouped.regions.table$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])


## Read in prior results
### Roberts et al 2021
jones_2012 <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Jones-et-al-2012-CSS-02-vDUKE.bed"))

jones_2012 <- jones_2012 %>%
  mutate(chr = factor(as.character(as.roman(gsub("chr", "", V1))), levels = levels(CSS.long$chr)),
        start = V2, end = V3, CSS_value = V4, split_num = V5) %>%
        select(-V1, -V2, -V3,-V4, -V5)

jones_2012$start.cum <- jones_2012$start+(chr$Cum.Seq.length[match(jones_2012$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])
jones_2012$end.cum <- jones_2012$end+(chr$Cum.Seq.length[match(jones_2012$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])

### Jones et al 2012
Roberts_2021 <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Prior_gasAcu-results/Roberts-et-al-2021-Specific-EcoPeaks-vDUKE.bed"))

Roberts_2021 <- Roberts_2021 %>%
  mutate(chr = factor(as.character(as.roman(gsub("chr", "", V1))), levels = levels(CSS.long$chr)),
        start = V2, end = V3) %>%
        select(-V1, -V2, -V3,-V4, -V5)

Roberts_2021$start.cum <- Roberts_2021$start+(chr$Cum.Seq.length[match(Roberts_2021$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])
Roberts_2021$end.cum <- Roberts_2021$end+(chr$Cum.Seq.length[match(Roberts_2021$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])

## Create ranges objects
Roberts_2021_ranges <- IRanges(Roberts_2021$start.cum, Roberts_2021$end.cum)
jones_2012_ranges <- IRanges(jones_2012$start.cum, jones_2012$end.cum)
top.grouped.regions.table_ranges <- IRanges(top.grouped.regions.table$start.cum, top.grouped.regions.table$end.cum)

## Calculate with regions of significance overlap with Jones and Roberts
top.grouped.regions.table$Overlap.Jones2012 <- FALSE
top.grouped.regions.table$Overlap.Jones2012[queryHits(findOverlaps(top.grouped.regions.table_ranges, jones_2012_ranges))] <- TRUE

top.grouped.regions.table$Overlap.Roberts <- FALSE
top.grouped.regions.table$Overlap.Roberts[queryHits(findOverlaps(top.grouped.regions.table_ranges, Roberts_2021_ranges))] <- TRUE

# Write out
top.grouped.regions.table %>%
  select(chr, start, end, mn.CSS, mn.CSS.sig, Overlap.Jones2012, Overlap.Roberts, genes) %>%
  write.table(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_top_regions_grouped.txt"), row.names = F, quote = F, sep = ",")

## Create cumulative position
chr$Sequence.name <- factor(chr$Sequence.name, levels = levels(CSS.long$chr))
chr <- chr[order(chr$Sequence.name),]
# Calculate cumulative length (starting at 0 )
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))
# Amened onto CSS data
CSS.long$start.cum <- CSS.long$start+(chr$Cum.Seq.length[match(CSS.long$chr, factor(chr$Sequence.name, levels = levels(CSS.long$chr)))]+chr$Cum.Seq.length[1])
top.regions.table$start.cum <- top.regions.table$start+(chr$Cum.Seq.length[match(top.regions.table$chr, factor(chr$Sequence.name, levels = levels(CSS.long$chr)))]+chr$Cum.Seq.length[1])

library(ggrepel)

# Plot
p.sig.windows <- ggplot(CSS.long) +
  geom_vline(xintercept = chr$Cum.Seq.length, col = "grey80") +
  geom_point(aes(start.cum, css, col = qval.0001<0.0001), size = 0.5) +
  # facet_grid(dropped~.) +
  geom_text_repel(data = top.regions.table[order(top.regions.table$mx.CSS, decreasing = T)[1:15],], 
            aes(x = start.cum, y = mx.CSS, label = contains.genes), nudge_y = 1, direction = "both", box.padding = 1) + #,hjust = -1, vjust = -1) +
  scale_color_manual(values=c("grey50", "deepskyblue"), name = "qvalue\nsignificance\n(<0.0001)") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Cum.Seq.length),20e6)),name = "Mbs", expand = c(0,0)) +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

# ggsave("test.png", p.sig.windows, width = 7.96*2, height = 7.96/2)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_annotated.png"), p.sig.windows, width = 7.96*2, height = 7.96/2)

p.sig.windows <- ggplot(CSS.long) +
  geom_vline(xintercept = chr$Cum.Seq.length, col = "grey80") +
  geom_point(aes(start.cum, css, col = qval.0001<0.0001), size = 0.5) +
  facet_grid(dropped~.) +
  geom_text_repel(data = top.regions.table[order(top.regions.table$mx.CSS, decreasing = T)[1:15],], 
            aes(x = start.cum, y = mn.CSS, label = contains.genes), size = 2, nudge_y = 1, nudge_x = 100, direction = "both", box.padding = 1) + #,hjust = -1, vjust = -1) +
  scale_color_manual(values=c("grey50", "deepskyblue"), name = "qvalue\nsignificance\n(<0.0001)") +
  theme_bw() +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Cum.Seq.length),20e6)),name = "Mbs", expand = c(0,0)) +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank(), legend.position = "bottom")

# ggsave("test.png", p.sig.windows, width = 7.96*2, height = 7.96*2)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_annotated_bydropPop.png"), p.sig.windows, width = 7.96*2, height = 7.96*2)


######################################
# # # # # DropPop CSS stats # # # # #
######################################

print("XXXXXXX PRINT STATS XXXXXXXXX")

CSS.long.all.sig.sum.stats <- CSS.long.all.sig %>% 
  group_by(dropped) %>%
  summarise(mn.CSS = mean(css, na.rm = T),
            md.CSS = median(css, na.rm = T),
            nm.wnd = sum(!is.na(css)),
            nm.sig.drop = sum(na.omit(qval.0001<0.0001)),
            nm.sig.all = sum(na.omit(all.sig.qvalue.0001))) %>%
  mutate(perc.sig.drop = (nm.sig.drop/nm.wnd)*100)

CSS.long.all.sig.sum.stats

top.regions.table %>%
  group_by() %>%
  summarise(nm.contig.wnd = length(chr),
            mn.length = mean(wnd.length, na.rm = T),
            mx.wnd.length = max(wnd.length),
            min.wnd.length = min(wnd.length),
            md.dist = median(bp.break[!is.infinite(bp.break)]),
            nm.within100kp = sum(bp.break[!is.infinite(bp.break)]<100000))


print("Sum of sig window legnths")
sum(top.regions.table$wnd.length)
(sum(top.regions.table$wnd.length)/sum(chr$Seq.length))*100

top.regions.table.chr <- top.regions.table %>%
  group_by(chr) %>%
  summarise(nm.contig.wnd = length(chr),
            nm.within100kp = sum(bp.break[!is.infinite(bp.break)]<100000),
            sum.length = sum(wnd.length, na.rm = T),) %>%
  mutate(nm.contig.within100kp = nm.contig.wnd-nm.within100kp) %>%
  select(-nm.within100kp)

top.regions.table.chr %>%
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_sig_window_bychr.csv"), row.names = F, quote = F, sep = ",")

top.regions.table.chr$chr[order(top.regions.table.chr$sum.length, decreasing = T)]
top.regions.table.chr$chr[order(top.regions.table.chr$nm.contig.within100kp, decreasing = T)]

# # # # # # # # # # # # # # # #
### Venn diagram of regions ###
# # # # # # # # # # # # # # # #

library(ggVennDiagram)

## Add cumulative length to CSS.wide
CSS.wide$start.cum <- CSS.wide$start+(chr$Cum.Seq.length[match(CSS.wide$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])
CSS.wide$end.cum <- CSS.wide$end+(chr$Cum.Seq.length[match(CSS.wide$chr, chr$Sequence.name)]+chr$Cum.Seq.length[1])

## Double check no duplicates
any(duplicated(CSS.wide$start.cum))

## Create ranges objects
Roberts_2021_ranges <- IRanges(Roberts_2021$start.cum, Roberts_2021$end.cum)
jones_2012_ranges <- IRanges(jones_2012$start.cum, jones_2012$end.cum)
CSS.wide_ranges <- IRanges(CSS.wide$start.cum, CSS.wide$end.cum)

## Calculate with regions of significance overlap with Jones and Roberts
CSS.wide$Overlap.Jones2012 <- FALSE
CSS.wide$Overlap.Jones2012[queryHits(findOverlaps(CSS.wide_ranges, jones_2012_ranges))] <- TRUE

CSS.wide$Overlap.Roberts <- FALSE
CSS.wide$Overlap.Roberts[queryHits(findOverlaps(CSS.wide_ranges, Roberts_2021_ranges))] <- TRUE

## Save CSS with overlaps
CSS.wide %>%
  write.table(file = paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_combine_with_overlaps.csv"), row.names = F, quote = F, sep = ",")

## Calc Venn Overlaps
venn.list <- list(This_study = CSS.wide$start.cum[CSS.wide$all.sig.qvalue.0001],
                  Jones_et_al_2012 = CSS.wide$start.cum[CSS.wide$Overlap.Jones2012],
                  Roberts_et_al_2021 = CSS.wide$start.cum[CSS.wide$Overlap.Roberts])


print("XXXXXXX PRINT STATS XXXXXXXX")
print("How many 2500 window were significant in Jones and this study")
length(CSS.wide$start.cum[CSS.wide$Overlap.Jones2012&!is.na(CSS.wide$all.sig.qvalue.0001)])
print("How many 2500 window were significant in Roberts and this study")
length(CSS.wide$start.cum[CSS.wide$Overlap.Roberts&!is.na(CSS.wide$all.sig.qvalue.0001)])

# Create plot
pVenn <- ggVennDiagram(venn.list) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")


## Venn Diagram for each individual dropped population
venn.list.lagoons <- list(nCLAC = CSS.wide$start.cum[CSS.wide$qval.sig.0001_nCLAC&!is.na(CSS.wide$qval.sig.0001_nCLAC)],
                  nLUIB = CSS.wide$start.cum[CSS.wide$qval.sig.0001_nLUIB&!is.na(CSS.wide$qval.sig.0001_nLUIB)],
                  nOBSE = CSS.wide$start.cum[CSS.wide$qval.sig.0001_nOBSE&!is.na(CSS.wide$qval.sig.0001_nOBSE)],
                  nDUIN = CSS.wide$start.cum[CSS.wide$qval.sig.0001_nDUIN&!is.na(CSS.wide$qval.sig.0001_nDUIN)],
                  Roberts_et_al_2021 = CSS.wide$start.cum[CSS.wide$Overlap.Roberts])

pVenn.lag <- ggVennDiagram(venn.list.lagoons) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")

# ggsave(paste0("test.png"), pVenn + pVenn.lag, width = 18, height = 9)
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_venn_studies.png"), pVenn + pVenn.lag, width = 19, height = 10)

## Venn Diagram for each individual dropped population
venn.list.HQ <- list(This_study = CSS.wide$start.cum[CSS.wide$all.sig.qvalue.0001&!is.na(CSS.wide$all.sig.qvalue.0001)],
                  Jones_et_al_2012 = CSS.wide$start.cum[CSS.wide$Overlap.Jones2012&!is.na(CSS.wide$all.sig.qvalue.0001)],
                  Roberts_et_al_2021 = CSS.wide$start.cum[CSS.wide$Overlap.Roberts&!is.na(CSS.wide$all.sig.qvalue.0001)])

pVenn.HQ <- ggVennDiagram(venn.list.HQ) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")

ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_venn_studies_HQ.png"), pVenn.HQ , width = 18, height = 9)

venn.list.lagoons.HQ.ndrop <- list(nCLAC = CSS.long$start.cum[CSS.long$dropped=="nCLAC"&CSS.long$qval.sig.0001&!is.na(CSS.long$qval.sig.0001)&CSS.long$css>=2],
                  nLUIB = CSS.long$start.cum[CSS.long$dropped=="nLUIB"&CSS.long$qval.sig.0001&!is.na(CSS.long$qval.sig.0001)&CSS.long$css>=2],
                  nOBSE = CSS.long$start.cum[CSS.long$dropped=="nOBSE"&CSS.long$qval.sig.0001&!is.na(CSS.long$qval.sig.0001)&CSS.long$css>=2],
                  nDUIN = CSS.long$start.cum[CSS.long$dropped=="nDUIN"&CSS.long$qval.sig.0001&!is.na(CSS.long$qval.sig.0001)&CSS.long$css>=2])

pVenn.lag.HQ.ndrop <- ggVennDiagram(venn.list.lagoons.HQ.ndrop) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")

ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_venn_studies_DropPop_HQ.png"), pVenn.lag.HQ.ndrop, width = 18, height = 9)

## Recalculate with regions of significance overlap with Jones and Roberts but removing large putative inversions
# Subset windows to outside of inversions
# filter split Iranges object to those smaller than 5Mb

p.hist <- ggplot(data.frame(width=do.call("c", width(top.regions.table.redueced)))) +
  geom_histogram(aes(x = width), bins = 50)

# Subset to those greater than 100kb
top.regions.table.redueced.Inv <- top.regions.table.redueced[(width(top.regions.table.redueced)>10000)]

## Filter out those regions from CSS.long
top.regions.table.redueced.Inv.hits <- findOverlaps(split(IRanges(CSS.wide$start, CSS.wide$end), CSS.wide$chr), top.regions.table.redueced.Inv)
# Remove these from CSS.wide
CSS.wide.noInv <- CSS.wide[-unique(queryHits(top.regions.table.redueced.Inv.hits)),]

# Recalculate Venn Overlaps
## Calc Venn Overlaps
venn.list.noInv <- list(This_study = CSS.wide.noInv$start.cum[CSS.wide.noInv$all.sig.qvalue.0001&!is.na(CSS.wide.noInv$all.sig.qvalue.0001)],
                  Jones_et_al_2012 = CSS.wide.noInv$start.cum[CSS.wide.noInv$Overlap.Jones2012&!is.na(CSS.wide.noInv$all.sig.qvalue.0001)],
                  Roberts_et_al_2021 = CSS.wide.noInv$start.cum[CSS.wide.noInv$Overlap.Roberts&!is.na(CSS.wide.noInv$all.sig.qvalue.0001)])

# Create plot
pVenn.noInv <- ggVennDiagram(venn.list.noInv) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")

# Save plot
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_venn_studies_noInv.png"), pVenn.noInv, width = 9, height = 9)


# Recalculate Venn Overlaps
## Calc Venn Overlaps
venn.list.minCSS2 <- list(This_study = CSS.wide$start.cum[CSS.wide$all.sig.qvalue.0001&!is.na(CSS.wide$all.sig.qvalue.0001)&CSS.wide$mn.CSS>=2],
                  Jones_et_al_2012 = CSS.wide$start.cum[CSS.wide$Overlap.Jones2012&!is.na(CSS.wide$all.sig.qvalue.0001)],
                  Roberts_et_al_2021 = CSS.wide$start.cum[CSS.wide$Overlap.Roberts&!is.na(CSS.wide$all.sig.qvalue.0001)])

# Create plot
pVenn.minCSS2 <- ggVennDiagram(venn.list.minCSS2) + scale_fill_gradient(low="grey90",high = "red") + theme(legend.position = "none")

# Save plot
ggsave(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_venn_studies_minCSS2.png"), pVenn.noInv, width = 9, height = 9)

