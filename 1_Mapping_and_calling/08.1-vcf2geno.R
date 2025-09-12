library(vcfR)
library(LEA)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
# vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.vcf.gz"
# vcf.ver <- "GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20"
## VCF version
vcf.ver <- args[2]

# Remove file extension
SNP.library.name <- basename(gsub(".vcf.gz", "", vcf.file))

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
  vcf.dir <- paste0("/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/", vcf.ver)
}

## Read in vcf file
vcf.SNPs <- read.vcfR(vcf.file, verbose = F)
vcf.SNPs

## Remove multiallelic snps and snps that are nolonger polymorphic
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs,na.omit = T),]

# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 
# Get sample names
samples <- colnames(vcf.SNPs@gt)[-1]

print("Converting vcf to geno")
geno.mat <- extract.gt(vcf.SNPs, return.alleles = F)
# Are any genotype calls not valid?
any(!unique(c(geno.mat))%in%c("0/0","0/1", "1/1", NA))

# Convert allels to numbers
geno.mat[geno.mat=="1/1"] <- 2
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="0/0"] <- 0

# Check none of the SNPs are entirely heterozgous and remove them if they are
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}

# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9
# Convert to numeric for storage
geno.mat.R <- apply(t(geno.mat), MARGIN = 1, as.numeric)

print("Writing out geno file.")
write.geno(geno.mat.R, paste0(vcf.dir, "/",SNP.library.name,".geno"))
#Write out sample names in order of geno
writeLines(samples, paste0(vcf.dir, "/",SNP.library.name,".geno.samples"))

## Create .phy file for RAxML
gt <- extract.gt(vcf.SNPs, return.alleles = T)

#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes,
gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- "R"
gt[gt=="G/A"] <- "R"
gt[gt=="C/T"] <- "Y"
gt[gt=="T/C"] <- "Y"
gt[gt=="A/C"] <- "M"
gt[gt=="C/A"] <- "M"
gt[gt=="G/T"] <- "K"
gt[gt=="T/G"] <- "K"
gt[gt=="C/G"] <- "S"
gt[gt=="G/C"] <- "S"
gt[gt=="A/T"] <- "W"
gt[gt=="T/A"] <- "W"
gt[gt=="."] <- "-"

# Check none of the SNPs are entirely heterozgous and remove them if they are
no.longer.poly <- apply(gt, MARGIN = 1, function(x) length(unique(x[x!="-"]))>1)
gt <- gt[no.longer.poly,]

## Write out file
ape::write.dna(t(gt), file = paste0(vcf.dir, "/",SNP.library.name,".phy"), format ="interleaved")

## Write out phy file with abigous bases removed
gt[!gt%in%c("A","T","C","G","-")] <- "-"

no.longer.poly <- apply(gt, MARGIN = 1, function(x) length(unique(x[x!="-"]))>1)
gt <- gt[no.longer.poly,]

ape::write.dna(t(gt), file = paste0(vcf.dir, "/",SNP.library.name,"nAbig.phy"), format ="interleaved")

