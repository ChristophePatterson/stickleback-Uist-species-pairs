library(vcfR)
library(LEA)
library(data.table)
library(tidyverse)

# Read in vcf
args <- commandArgs(trailingOnly=T)
# vcf.file <- args[1]
vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_DUIN_chrI_inv.vcf.gz"
vcf.SNPs <- read.vcfR(vcf.file)
outdir <- dirname(vcf.file)

## Get coding region of atp1a1a
# Read vcf
# vcf.atp1a1a.file <- args[2]
vcf.atp1a1a.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_DUIN_chrI_inv_SNPs_atp1a1a.vcf.gz"
vcf.atp1a1a <- read.vcfR(vcf.atp1a1a.file)

# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[, c("FORMAT", sort(colnames(vcf.SNPs@gt)[-1]))] 

## Get list of chromosomes
chr <- unique(vcf.SNPs@fix[,"CHROM"])

# Set window size
wndsize <- 5000
wndslid <- 2500

# Create windows for each chromosome
## Blank list
chr.df <- as.data.table(vcf.SNPs@fix)[, .(POS = as.numeric(POS)), by = CHROM]
sldwindows.df <- chr.df[, {
  end <- seq(min(POS), max(POS), by = wndslid)
  start <- end - wndsize + 1
  .(start = start, end = end)
}, by = CHROM]

# Calculate number of snps within each window
sldwindows.df$nsnps <- apply(sldwindows.df, MARGIN = 1, function(x) sum(between(as.numeric(chr.df$POS), as.numeric(x["start"]), as.numeric(x["end"]))))

geno.mat.full <- extract.gt(vcf.SNPs, element = "GT")

print("Starting PCA and MDS run")

mds.comp <- map2_dfr(
    sldwindows.df$CHROM, seq_len(nrow(sldwindows.df)),
    ~{
      chr <- .x
      idx <- .y
      # chr <- "CM102076.1"
      # idx <- 20
      wnd <- sldwindows.df[idx]
      wndname <- paste0(chr,"-", wnd$start,"-",wnd$end)


      snp.idx <- which(vcf.SNPs@fix[, "CHROM"] == chr &
                       between(as.numeric(vcf.SNPs@fix[, "POS"]), wnd$start, wnd$end))
      print(paste(wndname, "nsnps",length(snp.idx)))
      if(length(snp.idx) < 2) { return(NULL) ; print("Not enough SNPs or samples")}
      ## Set geno to be 0,1, or 9
      geno.mat <- geno.mat.full[snp.idx, , drop = FALSE]
      geno.mat <- matrix(as.integer(factor(geno.mat, levels = c("0/0", "0/1", "1/1")))-1, 
                     nrow = nrow(geno.mat), dimnames = dimnames(geno.mat))
      geno.mat <- apply(geno.mat, 2, as.integer)
      
      # Filter rows/cols with excessive NA, constant values, etc.
      is.samp.all.missing <- apply(geno.mat, MARGIN = 2, function(x) sum(is.na(x))>nrow(geno.mat)*0.5)
      if(any(is.samp.all.missing)){geno.mat <- geno.mat[,-is.samp.all.missing]}

      ##is.bad.snp <- apply(geno.mat, 1, function(x) {
      ##      vals <- unique(na.omit(x))
      ##      length(vals) <= 1 || all(vals == 1)  # All het or constant
      ##})
      ## if(any(is.bad.snp)){geno.mat <- geno.mat[-which(is.bad.snp), ]}
      print(paste("Stage 1", nrow(geno.mat), ncol(geno.mat)))
      #Check is matrix still contains data
      if(ncol(geno.mat)<=10|nrow(geno.mat)<=10){ return(NULL) ; print("Not enough SNPs or samples")}
      # Run MDS / PCA on cleaned geno.sub
      # print("Stage 2")
      ## Calculate distance matrix
      dc <- dist(t(geno.mat))
      dc.marine.samp <- as.matrix(dc)["Uist22617",]
      dc.marine.samp <- dc.marine.samp[colnames(geno.mat)]
      # Calculate number of called bases per sample
      base.calls <- apply(geno.mat, MARGIN = 2, function(x) sum(!is.na(x)))[colnames(geno.mat)]

      # Check is dist matrix contains any NA
      if(any(is.na(dc))) {return(NULL) ; print("mds-failed")}
      mds <- cmdscale(dc, k = 2)

      # Return data.frame with PCA1, PCA2, MDS1, MDS2, sample names, window name
      return(data.frame(sample = colnames(geno.mat), chr = chr, start = wnd$start, end = wnd$end, base.calls = base.calls, nsamps = ncol(geno.mat),
                        MDS1 = mds[,1], MDS2 = mds[,2], dist.marine = dc.marine.samp))

    }
  )


samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(mds.comp, samples,  by.y = "individual", by.x="sample", all.x = T)

# Read in atp1a1a
atp1a1a.bed <- read.table(paste0(outdir, "/atp1a1a_CDS.bed")) %>%
    rename(chr = V1, start = V2, end = V3)

p <- ggplot(samples_data) +
    geom_line(aes(start, dist.marine/base.calls, col = Population, group = sample)) +
    geom_vline(xintercept = c(26836909, 26867066)) +
    annotate("text", x = (26836909+((26867066-26836909)/2)), y = 0.0, label = "atp1a1a", color = "black", size = 4, vjust = 1) +
    geom_segment(data = atp1a1a.bed, aes(x = start, xend = end, y = -0.0001), linewidth = 10) +
    theme_bw() +
    coord_cartesian(clip = "off")  # allow drawing in the margin area

## ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"), p, width = 20)
## ggsave(paste0("test.png"), p, width = 20)

# Extract atp1a1a genotype
geno.mat.atp1a1a <- extract.gt(vcf.atp1a1a, element = "GT")
# Convert to geno
geno.mat.atp1a1a <- matrix(as.integer(factor(geno.mat.atp1a1a , levels = c("0/0", "0/1", "1/1")))-1, 
                     nrow = nrow(geno.mat.atp1a1a ), dimnames = dimnames(geno.mat.atp1a1a ))
geno.mat.atp1a1a  <- apply(geno.mat.atp1a1a , 2, as.integer)

# Calc dist
dc.atp1a1a <- dist(t(geno.mat.atp1a1a))
dc.atp1a1a.marine.samp <- as.matrix(dc.atp1a1a)["Uist22617",]
dc.atp1a1a.marine.samp <- dc.atp1a1a.marine.samp[colnames(geno.mat.atp1a1a)]
# Calculate number of called bases per sample
base.calls.atp1a1a <- apply(geno.mat.atp1a1a, MARGIN = 2, function(x) sum(!is.na(x)))[colnames(geno.mat.atp1a1a)]

# Run mds
mds.atp1a1a <- cmdscale(dc.atp1a1a, k = 2)

# Combine
mds.atp1a1a.df <- data.frame(sample = colnames(geno.mat.atp1a1a), base.calls = base.calls.atp1a1a, nsamps = ncol(geno.mat.atp1a1a),
                        MDS1 = mds.atp1a1a[,1], MDS2 = mds.atp1a1a[,2], dist.marine = dc.atp1a1a.marine.samp) %>%
                         merge(samples,  by.y = "individual", by.x="sample", all.x = T)

#  Calc nj 
library(ape)
library(ggtree)
library(patchwork)
nj.data <- nj(dc.atp1a1a)

# Create tree plot
plot.tree <- ggtree(nj.data, aes(color = Population), layout = "daylight")
#  Combine with sample data
plot.tree <- plot.tree %<+% mds.atp1a1a.df
## Custom tip colours
plot.tree <- plot.tree + geom_tippoint(aes(fill = Population, shape = Ecotype), stroke = 1, size=3) 

ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"), p + plot.tree + plot_layout(widths = c(3,1)), width = 20)
ggsave(paste0("test.png"), p + plot.tree + plot_layout(widths = c(3,1)), width = 20)
