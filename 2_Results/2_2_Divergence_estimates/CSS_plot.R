library(ggplot2)
## devtools::install_github("jdstorey/qvalue")
library(qvalue)


args<-commandArgs(trailingOnly=T)
CSS.dir <- args[1]
CSS.run <- args[2]

chr <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv",
                  header = T, sep = "\t")

setwd(CSS.dir)

CSS <- read.table(paste0(CSS.dir,"/", CSS.run),header = T)

chr <- chr[!chr$Chromosome.name%in%c("Un","MT"),]

chr$Chromosome.name <- factor(chr$Chromosome.name, levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

CSS$chr <- factor(chr$Chromosome.name[match(CSS$chr, chr$RefSeq.seq.accession)], levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

CSS.q <- qvalue(1-CSS$pval, fdr.level = 0.0001)

summary(CSS.q)
## Filter out CSS
CSS$qval.calc <- CSS.q$qvalues
CSS$qval.sig <- CSS.q$significant
dim(CSS)

p <- ggplot(CSS) +
  geom_point(aes(start, css, col = qval.sig)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  scale_color_manual(values=c("grey50", "deepskyblue")) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines')) +
  ggtitle(gsub(".txt","",CSS.run))

ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),".pdf"), p, height = 5, width = 15)
ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),".png"), p, height = 5, width = 15)

p <- ggplot(CSS) +
  geom_point(aes(start, css, col = qval.sig)) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  scale_color_manual(values=c("grey50", "deepskyblue")) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6), breaks = c(seq(0, max(chr$Seq.length),10e6)),name = "Mbs") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'))

ggsave(paste0(CSS.dir,"/",gsub(".txt","",CSS.run),"_vert.pdf"), p, height = 30, width = 15)
ggsave(paste0(CSS.dir,"/", gsub(".txt","",CSS.run),"_vert.png"), p, height = 30, width = 15)


