
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
library(data.table)

#Not currently installed
#library(poppr)
#library(ggnewscale)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
# vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/pca/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX/NC_053212.1/stickleback.NC_053212.1.vcf.gz"
vcf.ver <- args[2]
# vcf.ver <- "ploidy_aware_HWEPops_MQ10_BQ20"
wndsize <- as.numeric(args[3])
# wndsize <- 25000
wndslid <- as.numeric(args[4])
# wndslid <- 5000
run_analysis <- args[5]

## vcf.file <- "stickleback_SNPs.rand10000.vcf.gz"
# Remove file extension
SNP.library.name <- basename(gsub(".vcf.gz", "", vcf.file))
plot.dir <- gsub(paste0(SNP.library.name,".vcf.gz"), "", vcf.file)

## Create directory is not already
dir.create(plot.dir)

# Set workign directory
setwd(plot.dir)

if(run_analysis){
  ## Read in vcf file
  vcf.SNPs <- read.vcfR(vcf.file, verbose = F)
  # Make vcf be in alphabetical order
  vcf.SNPs <- vcf.SNPs[, c("FORMAT", sort(colnames(vcf.SNPs@gt)[-1]))] 
  
  # Get an read sample information
  sample_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)

  ######################################
  ##### PCA for all samples ####
  ######################################
  ## Remove multiallelic snps and snps that are nolonger polymorphic
  vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
  vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs, na.omit = T),]

  samples <- data.frame(samples = colnames(vcf.SNPs@gt)[-1])
  ## Get list of chromosomes
  chr <- unique(vcf.SNPs@fix[,"CHROM"])

  # Create windows for each chromosome
  ## Blank list

  chr.df <- as.data.table(vcf.SNPs@fix)[, .(POS = as.numeric(POS)), by = CHROM]
  sldwindows.df <- chr.df[, {
    end <- seq(wndsize, max(POS), by = wndslid)
    start <- end - wndsize + 1
    .(start = start, end = end)
  }, by = CHROM]

  # Calculate number of snps within each window
  sldwindows.df$nsnps <- apply(sldwindows.df, MARGIN = 1, function(x) sum(between(as.numeric(chr.df$POS), x["start"], x["end"])))

  geno.mat.full <- extract.gt(vcf.SNPs, element = "GT")

  print("Starting PCA and MDS run")
  pca.comp <- map2_dfr(
    sldwindows.df$CHROM, seq_len(nrow(sldwindows.df)),
    ~{
      chr <- .x
      idx <- .y
      wnd <- sldwindows.df[idx]
      wndname <- paste0(chr,"-", wnd$start,"-",wnd$end)

      snp.idx <- which(vcf.SNPs@fix[, "CHROM"] == chr &
                       between(as.numeric(vcf.SNPs@fix[, "POS"]), wnd$start, wnd$end))

      if(length(snp.idx) < 2) return(NULL)

      geno.mat <- geno.mat.full[snp.idx, , drop = FALSE]
      geno.mat <- matrix(as.integer(factor(geno.mat, levels = c("0/0", "0/1", "1/1")))-1, 
                     nrow = nrow(geno.mat), dimnames = dimnames(geno.mat))
      geno.mat <- apply(geno.mat, 2, as.integer)

      # Filter rows/cols with excessive NA, constant values, etc.
      is.samp.all.missing <- apply(geno.mat, MARGIN = 1, function(x) sum(is.na(x))>ncol(geno.mat)*0.5)
      if(any(is.samp.all.missing)){geno.mat <- geno.mat[-which(is.samp.all.missing),]}

      is.bad.snp <- apply(geno.mat, 2, function(x) {
            vals <- unique(na.omit(x))
            length(vals) <= 1 || all(vals == 1)  # All het or constant
      })
      geno.mat <- geno.mat[, !is.bad.snp, drop = FALSE]

      #Check is matrix still contains data
      if(ncol(geno.mat)<=10|nrow(geno.mat)<=10) return(NULL)
      # Run MDS / PCA on cleaned geno.sub

      ## Calculate distance matrix
      dc <- dist(t(geno.mat))
      # Check is dist matrix contains any NA
      if(any(is.na(dc))) return(NULL)
      mds <- cmdscale(dc, k = 2)

      # Make missing SNPs equal to "9"
      geno.mat[is.na(geno.mat)] <- 9

      geno.df <- data.frame(geno.mat)
      dim(geno.df)

      print("Writing out geno file.")
      write.table(x = geno.df, file = paste0(plot.dir, SNP.library.name, wndname,".geno"),
                      col.names = F, row.names = F, quote = F, sep = "")

      #Conduct PCA
      print("Conducting PCA")
      geno2lfmm(paste0(plot.dir, SNP.library.name, wndname,".geno"), 
                    paste0(plot.dir, SNP.library.name, wndname,".lfmm"), force = TRUE)
      #PCA
      print("Reading back in PCA")
      pc <- pca(paste0(plot.dir, SNP.library.name, wndname,".lfmm"), scale = TRUE)

      # Return data.frame with PCA1, PCA2, MDS1, MDS2, sample names, window name
      return(data.frame(sample = colnames(geno.mat), chr = chr, start = wnd$start, end = wnd$end, nsnps = nrow(geno.mat), nsamps = ncol(geno.mat),
                        PCA1 = pc$projections[,1], PCA2 = pc$projections[,2], MDS1 = mds[,1], MDS2 = mds[,1]))


    }
  )
  write.csv(pca.comp, paste0(plot.dir, SNP.library.name, "_sliding-window_pca_wndsize", format(wndsize,scientific = F),"_wndslid", format(wndslid, scientific = F),".txt"),
            row.names = F)

  }

print("PCA and MDS analysis completed")
