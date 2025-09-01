## devtools::install_github("jdstorey/qvalue")
library(qvalue)
library(tidyverse)


args<-commandArgs(trailingOnly=T)
CSS.dir <- args[1]
#CSS.dir <-  "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/stickleback.wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05"
CSS.run <- args[2]
#CSS.run <- "stickleback.wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05.CSSm.10000perm.txt"

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

setwd(CSS.dir)

CSS <- read.table(paste0(CSS.dir,"/", CSS.run),header = T)

chr <- chr[!chr$Chromosome.name%in%c("Un","MT"),]

chr <- chr[chr$GenBank.seq.accession%in%unique(CSS$chr),]
chr$Sequence.name <- gsub("chr", "", chr$Sequence.name)
# chr <- chr[order(chr$Seq.length, decreasing = T),]
chr$bi.col <- rep(1:2, length.out = dim(chr)[1])

CSS$chr <- chr$Sequence.name[match(CSS$chr, chr$GenBank.seq.accession)]
CSS$bi.col <- chr$bi.col[match(CSS$chr, chr$GenBank.seq.accession)]

CSS$chr <- factor(CSS$chr, levels = chr$Sequence.name)

p.hist <- ggplot(CSS) +
  geom_histogram(aes(nsnps), bins = 100)

ggsave("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/test.png", p.hist)

## Remove low quality snps
CSS.HQ <- CSS[CSS$nsnps>=10&CSS$nsnps<=200,]
## Calculate q values
CSS.q.05 <- qvalue(1-CSS.HQ$pval, fdr.level = 0.05)
CSS.q.001 <- qvalue(1-CSS.HQ$pval, fdr.level = 0.001)
CSS.q.0001 <- qvalue(1-CSS.HQ$pval, fdr.level = 0.0001)
# Add to dataset
CSS.HQ$qval.calc.05 <- CSS.q.05$qvalues
CSS.HQ$qval.sig.05 <- CSS.q.05$significant
CSS.HQ$qval.calc.001 <- CSS.q.001$qvalues
CSS.HQ$qval.sig.001 <- CSS.q.001$significant
CSS.HQ$qval.calc.0001 <- CSS.q.0001$qvalues
CSS.HQ$qval.sig.0001 <- CSS.q.0001$significant

## Assign unique ID to each run of signifcant results (split by chromosome)
CSS.HQ <- CSS.HQ %>%
  mutate(temp = replace(qval.sig.05, is.na(qval.sig.05), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  mutate(goal.05 = ifelse(temp, paste0(first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

CSS.HQ <- CSS.HQ %>%
  mutate(temp = replace(qval.sig.001, is.na(qval.sig.001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  mutate(goal.001 = ifelse(temp, paste0(first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

CSS.HQ <- CSS.HQ %>%
  mutate(temp = replace(qval.sig.0001, is.na(qval.sig.0001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  mutate(goal.0001 = ifelse(temp, paste0(first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

## How many windows are there?
nrow(CSS.HQ)
## Range of SNPs
summary(CSS.HQ$nsnps)
mean(CSS.HQ$nsnps)
sd(CSS.HQ$nsnps)
# Mean CSS score
summary(CSS.HQ$css)
# Percentage were significant
(sum(CSS.HQ$qval.sig.05)/nrow(CSS.HQ))*100
(sum(CSS.HQ$qval.sig.001)/nrow(CSS.HQ))*100
(sum(CSS.HQ$qval.sig.0001)/nrow(CSS.HQ))*100

## Split into X number of regions
length(unique(CSS.HQ$goal.0001))
length(unique(CSS.HQ$goal.0001[CSS.HQ$css>=2]))
# Length distribution of signigicant windows
sort(table(CSS.HQ$goal.0001), decreasing = T)*2500
(sort(table(CSS.HQ$goal.0001), decreasing = T)*2500)[1:10]
summary(as.numeric(sort(table(CSS.HQ$goal.0001), decreasing = T)*2500))

## Top regions
CSS.HQ$goal.0001[CSS.HQ$css>=2]
top.regions <- sort(table(CSS.HQ$goal.0001[CSS.HQ$css>=2]), decreasing = T)*500
top.regions <- top.regions[top.regions>=(2500)]
length(top.regions)
CSS.HQ$top.regions <- CSS.HQ$goal.0001 %in% names(top.regions)
unique(CSS.HQ$chr[CSS.HQ$top.regions])

## Table
top.regions.table <- CSS.HQ %>%
  filter(top.regions) %>%
  group_by(goal.0001) %>%
  summarise(chr = first(chr),
            start = min(start),
            end = max(end),
            mn.CSS = mean(css),
            wnd.length = 1+max(end)-min(start)) %>%
  ungroup() %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(bp.break = c(Inf, diff(start)))

## Read in DUKE mapped bed annotations
DUKE.bed <- tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Duke_GAcu_1_CDS.gtf", header = F)) %>%
  rename(chr = V1, anon.type =  V2, start = V4, end = V5, type = V3) %>%
  mutate(chr = as.factor(as.character(as.roman(as.numeric(gsub("chr", "", chr))))))

DUKE.bed <- tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/stickleback_DUKE_ensembl_genes.bed", header = F)) %>%
  rename(chr = V1, start = V2, end = V3, type = V4, ID = V5) %>%
  mutate(chr = as.factor(as.character(as.roman(as.numeric(gsub("chr", "", chr))))))

DUKE.CSS.match <- DUKE.bed %>% 
  filter(chr %in% top.regions.table$chr) %>%
  filter(type == "gene") %>%
  rowwise() %>%
  filter(any(chr == top.regions.table$chr & start >= top.regions.table$start & end <= top.regions.table$end )) 

DUKE.CSS.match <- DUKE.CSS.match %>% 
  mutate(gene.name = gsub("Name=", "", str_split_i(ID, ";", 2)))

# Finds all genes located with and aroud the identified region
top.regions.table <- top.regions.table %>%
  rowwise() %>%
  mutate(contains.genes = paste0(DUKE.CSS.match$gene.name[DUKE.CSS.match$chr==chr&(between(DUKE.CSS.match$start, start, end)|between(DUKE.CSS.match$end, start, end))], collapse = ","),
         contains.genes.in10Kbps = paste0(DUKE.CSS.match$gene.name[DUKE.CSS.match$chr==chr&(between(DUKE.CSS.match$start, start-10000, end+10000)|between(DUKE.CSS.match$end, start-10000, end+10000))], collapse = ","),
         contains.genes.in100Kbps = paste0(DUKE.CSS.match$gene.name[DUKE.CSS.match$chr==chr&(between(DUKE.CSS.match$start, start-100000, end+100000)|between(DUKE.CSS.match$end, start-100000, end+100000))], collapse = ",")) %>%
  mutate(contains.genes.in100Kbps = gsub(pattern = contains.genes.in10Kbps, replacement = "", x = contains.genes.in100Kbps), # removes genes already included first gene contains column
         contains.genes.in10Kbps = gsub(pattern = contains.genes, replacement = "", x = contains.genes.in10Kbps))

# Write out file
top.regions.table %>% 
  select(-goal.0001, -bp.break) %>%
  mutate(mn.CSS = round(mn.CSS, 3)) %>%
  write.table(file = paste0(CSS.dir,"/",gsub(".txt","",CSS.run),"_top_regions.txt"), row.names = F, quote = F)

p <- ggplot(CSS.HQ) +
  geom_point(aes(start, css, col = qval.sig.0001)) +
  geom_point(data = CSS.HQ[CSS.HQ$top.regions,], aes(start, css), col = "firebrick") +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_color_manual(values=c("grey50", "deepskyblue")) +
  geom_errorbar(data = CSS.HQ[1,],
    aes(xmin = 2e6, xmax = 22e6, y = max(CSS.HQ$css)*0.8), linewidth = 1, width = 0.4)+
  geom_text(data = CSS.HQ[1,],
    aes(x = ((22e6-2e6)/2)+2e6, y = max(CSS.HQ$css)*0.82), label = "20Mbp", size = 3) +
  theme_classic() +
  theme(legend.position = "none",panel.spacing = unit(0,'lines')) +
  ggtitle(gsub(".txt","",CSS.run)) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  scale_y_continuous(sec.axis = sec_axis(~., labels = NULL, breaks = 0)) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        axis.line = element_line(), strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        axis.title.x = element_blank())

ggsave(paste0("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/test.png"), p, height = 5, width = 15)
ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),".pdf"), p, height = 5, width = 15)
ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),".png"), p, height = 5, width = 15)

p <- ggplot(CSS.HQ) +
  geom_point(aes(start, css, col = qval.sig.0001)) +
    geom_point(data = CSS.HQ[CSS.HQ$top.regions,], aes(start, css), col = "firebrick") +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  scale_color_manual(values=c("grey50", "deepskyblue")) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'))


print(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),"_vert.pdf"))
ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),"_vert.pdf"), p, height = 30, width = 15)
ggsave(paste0(CSS.dir,"/", gsub(".txt","",CSS.run),"_vert.png"), p, height = 30, width = 15)


