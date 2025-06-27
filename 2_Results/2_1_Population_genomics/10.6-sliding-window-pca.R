
# CalcuLation of LEA on stickleback popuLations
library(patchwork)
library(ggplot2)
library(ape)
library(vcfR)
library(tidyverse)
# BiocManager::install(version = '3.20')
# BiocManager::install("LEA")
library(LEA)
library(adegenet)
library(ggrepel)
library(scatterpie)

#Not currently installed
#library(poppr)
#library(ggnewscale)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
vcf.ver <- args[2]
wndsize <- as.numeric(args[3])
wndslid <- as.numeric(args[4])
run_analysis <- args[5]
## vcf.file <- "stickleback_SNPs.rand10000.vcf.gz"
# Remove file extension
SNP.library.name <- basename(gsub(".vcf.gz", "", vcf.file))

# Test whether working on HPC or laptop and set working directory accordingly
# HPC test
dir.path <-"/gpfs01/home/mbzcp2/data/sticklebacks/"
plot.dir <- paste0("/gpfs01/home/mbzcp2/data/sticklebacks/results/", vcf.ver, "/sliding-window/pca/", SNP.library.name, "/")
print(plot.dir)
# Set workign directory
setwd(plot.dir)

## Create directory is not already
dir.create(plot.dir)
## Read in vcf file
## Read in vcf file

if(run_analysis){
  vcf.SNPs <- read.vcfR(vcf.file, verbose = F)
  # Make vcf be in alphabetical order
  vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 

  ## Remove X and Y chromosomes
  sex_chr <- c("NC_053230.1","NC_053233.1")
  vcf.SNPs <- vcf.SNPs[!vcf.SNPs@fix[,1]%in%sex_chr]

  # Get an read sample information
  sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)

  ######################################
  ##### PCA for all samples ####
  ######################################
  ## Remove multiallelic snps and snps that are nolonger polymorphic
  vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
  vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs, na.omit = T),]

  ## Get list of chromosomes
  chr <- unique(vcf.SNPs@fix[,"CHROM"])

  # Create windows for each chromosome
  ## Blank list
  sldwindows <- list()
  ## Loop for each chromsome
  for(i in 1:length(chr)){
    # Get chromosome i
    chri <- chr[i]
    ## Get max position
    chri.max <- max(as.numeric(vcf.SNPs@fix[,"POS"])[vcf.SNPs@fix[,"CHROM"]==chri])
    # Get end locations
    end <- seq(wndsize, chri.max, by =  wndslid)
    # Get start locations
    start <- end-wndsize+1
    sldwindows[[i]] <- data.frame(chr = chri, start = start, end = end, nsnps = sum(vcf.SNPs@fix[,"CHROM"]==chri))
  }
  # Combine into dataframe
  sldwindows.df <- do.call("rbind", sldwindows)
  sldwindows.df[1,]
  dim(sldwindows.df)

  samples <- data.frame(samples = colnames(vcf.SNPs@gt)[-1])

  i <- 1
  pca.comp <- list()

  for(i in 1:length(sldwindows.df$chr)){
    wndname <- paste0(sldwindows.df$chr[i], "-",  format(sldwindows.df$start[i], scientific = F), "_", format(sldwindows.df$end[i], scientific = F))
    SNP.sub.tmp <- vcf.SNPs@fix[,"CHROM"]==sldwindows.df$chr[i]&between(as.numeric(vcf.SNPs@fix[,"POS"]), sldwindows.df$start[i], sldwindows.df$end[i])
    vcf.SNPs.window <- vcf.SNPs[SNP.sub.tmp]
    vcf.SNPs.window
    if(sum(SNP.sub.tmp)>=2){

      ## Convert to geno object
      print("Converting vcf to geno")
      geno <- vcfR2loci(vcf.SNPs.window, return.alleles = F)
      geno.mat <- as.matrix(geno)
      geno.mat[geno.mat=="1/1"] <- 2
      geno.mat[geno.mat=="0/1"] <- 1
      geno.mat[geno.mat=="0/0"] <- 0

      # Remove low covarage samples
      is.samp.all.missing <- apply(geno.mat, MARGIN = 1, function(x) sum(is.na(x))>dim(geno.mat)[2]*0.5)
      if(any(is.samp.all.missing)){geno.mat <- geno.mat[-which(is.samp.all.missing),]}

      # Check none of the SNPs are entirely heterozgous and remove them if they are
      is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
      if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}

      ## Check if any snps is contant across all individuals
      is.constant.snp <- apply(geno.mat, MARGIN = 2, function(x) nchar(gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = ""))<=1)
      if(any(is.constant.snp)){geno.mat <- geno.mat[,-which(is.constant.snp)]}

      print(dim(geno.mat))
      # Check if there are are still samples and genotypes
      if(dim(geno.mat)[2]<=10|dim(geno.mat)[1]<=10){print(paste0("Not enough SNPs in", wndname))} else {
        # Check none of the SNPs are entirely heterozgous and remove them if they are
        dim(geno.mat)

        ## MDS
        dc <- dist(geno.mat)
        mds <- cmdscale(dc, k = 2)      

        # Make missing SNPs equal to "9"
        geno.mat[is.na(geno.mat)] <- 9

        geno.df <- data.frame(t(geno.mat))
        dim(geno.df)

        print("Writing out geno file.")
        write.table(x = geno.df, file = paste0(plot.dir, wndname,".geno"),
                    col.names = F, row.names = F, quote = F, sep = "")

        #Conduct PCA
        print("Conducting PCA")
        geno2lfmm(paste0(plot.dir, wndname,".geno"), 
                  paste0(plot.dir, wndname,".lfmm"), force = TRUE)
        #PCA
        print("Reading back in PCA")
        pc <- pca(paste0(plot.dir, wndname,".lfmm"), scale = TRUE)

        # Stores PCA data
        pca.comp[[wndname]] <- data.frame(samples = colnames(geno.df), windowname = wndname, 
                                          PCA1 = pc$projections[,1], PCA2 = pc$projections[,2],
                                          MDS1 = mds[,1], MDS2 = mds[,2])
        # pca.comp[[wndname]] <- merge(samples, pca.comp.tmp, by = "samples", all.x = T)
        if (file.exists(paste0(plot.dir, wndname,".lfmm"))) {
          #Delete file if it exists
          file.remove(paste0(plot.dir, wndname,".lfmm"))
          file.remove(paste0(plot.dir, wndname,".geno"))
          file.remove(paste0(plot.dir, wndname,".pcaProject"))
          unlink(paste0(plot.dir, wndname,".pca"), recursive = T)
        }
      }
    }

    print(paste(round(i/length(sldwindows.df$chr)*100), "% of way through analysis"))
}

pca.comp.df <- do.call("rbind", pca.comp)
pca.comp.df$chr.name <- stringr::str_split_i(pca.comp.df$windowname, pattern = "-", 1)
pca.comp.df$end <- stringr::str_split_i(pca.comp.df$windowname, pattern = "_", 3)

## Add in Population, Ecotype, and Waterbody to local PCA
pca.comp.df$Population <- sample_data$Population[match(pca.comp.df$samples, sample_data$individual)]
pca.comp.df$Ecotype <- sample_data$Ecotype[match(pca.comp.df$samples, sample_data$individual)]
pca.comp.df$Waterbody <- sample_data$Waterbody[match(pca.comp.df$samples, sample_data$individual)]

# Save output
write.csv(pca.comp.df, paste0(plot.dir, "sliding-window_pca_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".txt"),
          row.names = F)

}

pca.comp.df <- read.csv(paste0(plot.dir, "sliding-window_pca_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".txt"), header = T)

# Read in chromosome details
chr <- chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T))
# Remove smaller scafs and mito
chr <- chr[!chr$Chromosome.name%in%c("Un","MT"),]
# Order Chromsome name to be in order 
chr$Chromosome.name <- factor(chr$Chromosome.name, levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

# Order chr in order in pca
pca.comp.df$chr <- chr$Chromosome.name[match(pca.comp.df$chr.name, chr$RefSeq.seq.accession)]
pca.comp.df$chr <- factor(pca.comp.df$chr, levels = chr$Chromosome.name[order(chr$RefSeq.seq.accession)])

# Transform PCA so that the axis is also segregating populations in the same direction across all windows
## Copy over PCA and MDS data to new scaled columns
pca.comp.df$MDS1_scaled <- pca.comp.df$MDS1
pca.comp.df$MDS2_scaled <- pca.comp.df$MDS2
pca.comp.df$PCA1_scaled <- pca.comp.df$PCA1
pca.comp.df$PCA2_scaled <- pca.comp.df$PCA2
## invert so population direction is consisten
for(wndow in unique(pca.comp.df$windowname)){
  ## Calculate if the median score for an Ecotype is positive or negative
  Pop.sign <- sign(median(pca.comp.df$MDS1_scaled[pca.comp.df$windowname==wndow&pca.comp.df$Ecotype=="anad"], na.rm=T))
  if(Pop.sign==(-1)){
    # If negative invert the score for all samples for that window
    pca.comp.df$MDS1_scaled[pca.comp.df$windowname==wndow] <- -pca.comp.df$MDS1_scaled[pca.comp.df$windowname==wndow]
  }
  # Repeat for MDS2
  Pop.sign <- sign(median(pca.comp.df$MDS2_scaled[pca.comp.df$windowname==wndow&pca.comp.df$Ecotype=="anad"], na.rm=T))
  if(Pop.sign==(-1)){
    pca.comp.df$MDS2_scaled[pca.comp.df$windowname==wndow] <- -pca.comp.df$MDS2_scaled[pca.comp.df$windowname==wndow]
  }
  # Repeat for PCA1
  Pop.sign <- sign(median(pca.comp.df$PCA1_scaled[pca.comp.df$windowname==wndow&pca.comp.df$Ecotype=="anad"], na.rm=T))
  if(Pop.sign==(-1)){
    pca.comp.df$PCA1_scaled[pca.comp.df$windowname==wndow] <- -pca.comp.df$PCA1_scaled[pca.comp.df$windowname==wndow]
  }
  # Repeat for PCA2
  Pop.sign <- sign(median(pca.comp.df$PCA2_scaled[pca.comp.df$windowname==wndow&pca.comp.df$Ecotype=="anad"], na.rm=T))
  if(Pop.sign==(-1)){
    pca.comp.df$PCA2_scaled[pca.comp.df$windowname==wndow] <- -pca.comp.df$PCA2_scaled[pca.comp.df$windowname==wndow]
  }
}


## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))
# Save output
ggsave(paste0(plot.dir, "sliding-window_pca_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p/q, width = 40, height = 15)

## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = samples)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), PCA2_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = samples)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))
# Save output
ggsave(paste0(plot.dir, "sliding-window_pca_line_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p/q, width = 40, height = 15)


## Plot just chromsome IX and I
p <- ggplot(pca.comp.df[pca.comp.df$chr=="IX"|pca.comp.df$chr=="I",], aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_wrap(~chr) +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

ggsave(paste0(plot.dir, "sliding-window_pca_IX_I_size", format(wndsize,scientific = F),"_slid", format(wndslid, scientific = F),".png"), p, width = 15, height = 10)

# Break down by chromosome
p <- ggplot(pca.comp.df, aes(as.numeric(end), PCA1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

ggsave(paste0(plot.dir, "sliding-window_pca_wndsize_Popsplit", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p, width = 40, height = 20)

## PLot MDS1 along genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

ggsave(paste0(plot.dir, "sliding-window_mds_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p, width = 40, height = 20)

## Plot MDS1 along genome and split by waterbody
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_point() +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

ggsave(paste0(plot.dir, "sliding-window_mds_wndsize_Popsplit", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p, width = 40, height = 20)

## Plot MDS1 along genome and split by waterbody
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = samples)) +
  facet_grid(Waterbody~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

ggsave(paste0(plot.dir, "sliding-window_mds_line_wndsize_Popsplit", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p, width = 40, height = 20)

## Plot PCA 1 across the genome
p <- ggplot(pca.comp.df, aes(as.numeric(end), MDS1_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = samples)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))

## Plot PCA 2 across the genome
q <- ggplot(pca.comp.df, aes(as.numeric(end), MDS2_scaled, col = Population, shape = Ecotype)) +
  geom_line(aes(group = samples)) +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(legend.position = "top",panel.spacing = unit(0,'lines'),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                      # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0)),  # move left axis labels closer to axis 
        strip.background = element_rect(size = 0.5))
# Save output
ggsave(paste0(plot.dir, "sliding-window_mds_line_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".png"), p/q, width = 40, height = 15)
