#Ada LEA and PCR plots
# Combined PCA and admixture plot

# Run 'module load R-uoneasy/4.2.1-foss-2022a' 

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
library(poppr)


#Not currently installed
#library(ggnewscale)
#library(treedataverse)

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
# vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/ploidy_aware_HWEPops_MQ10_BQ20/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.vcf.gz"
vcf.ver <- args[2]
## vcf.ver <- "ploidy_aware_HWEPops_MQ10_BQ20"
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
}
## Create directory is not already
dir.create(plot.dir)
dir.create(paste0(plot.dir, "/LEA_PCA/"))
dir.create(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/"))

## Read in vcf file
vcf.SNPs <- read.vcfR(vcf.file, verbose = T)
# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 

## Remove X and Y chromosomes
sex_chr <- c("NC_053230.1","NC_053233.1")
vcf.SNPs <- vcf.SNPs[!vcf.SNPs@fix[,1]%in%sex_chr]

colnames(vcf.SNPs@gt)

## Calculated sequecning error rate from duplicated samples
## technical_dups <- vcf.SNPs[samples = c("Obsm_640", "Obsm_641")]
## 
## technical_dups_gt <- extract.gt(technical_dups)
## sum(as.numeric(na.omit(!technical_dups_gt[,1]==technical_dups_gt[,2])))/length(na.omit(technical_dups_gt[,1]==technical_dups_gt[,2]))*100
## 
## technical_dups_gt_errors <- na.omit(technical_dups_gt[technical_dups_gt[,1]!=technical_dups_gt[,2],])
## table(paste0(technical_dups_gt_errors[,1], "-", technical_dups_gt_errors[,2]))

## Remove one of the dupicalted samples
vcf.SNPs <- vcf.SNPs[samples = colnames(vcf.SNPs@gt)[colnames(vcf.SNPs@gt)!="Obsm_641"]]

colnames(vcf.SNPs@gt)
# Get an read sample information
samples_data <- data.frame(ID = colnames(vcf.SNPs@gt)[-1])
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)
samples_data <- samples_data[samples_data$ID%in%colnames(vcf.SNPs@gt),]
samples_data <- samples_data[match(samples_data$ID, (colnames(vcf.SNPs@gt)[-1])),]
print("Do any samples names not line up with (False is good)")
any(!samples_data$ID==(colnames(vcf.SNPs@gt)[-1]))
## Fill in missing data
samples_data$Ecotype[is.na(samples_data$Ecotype)] <- "Unknown"

### Missing genotype assesment
# Stats
mySampleStats <- apply(extract.gt(vcf.SNPs), MARGIN = 2, function(x){ sum(is.na(x)) })
mySampleStats <- data.frame(sample = names(mySampleStats), per.gt = mySampleStats/nrow(vcf.SNPs)*100)
dim(mySampleStats)
mySampleStats <- merge(mySampleStats, samples_data[,c(1,7,8,9,10,11,12,13,14)], by.x = "sample", by.y = "ID")

summary(mySampleStats)
p <- ggplot(mySampleStats) +
  geom_boxplot(aes(x = Population, y = per.gt, col = Population), outlier.colour = NA) +
  geom_jitter(aes(x = Population, y = per.gt, col = Population), height = 0, width = 0.2) 

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_missGT.png"), p, width = 10, height = 10)

mySampleStats$mn_cov <- apply(extract.gt(vcf.SNPs, element = "DP"), MARGIN = 2, function(x){ mean(as.numeric(x)) })
mySampleStats$sd_cov <- apply(extract.gt(vcf.SNPs, element = "DP"), MARGIN = 2, function(x){ sd(as.numeric(x)) })

p <- ggplot(mySampleStats) +
  geom_boxplot(aes(x = Population, y = mn_cov, col = Population), outlier.colour = NA) +
  geom_jitter(aes(x = Population, y = mn_cov, col = Population), height = 0, width = 0.2) 

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_sampleCov.png"), p, width = 10, height = 10)

mySNPStats <- apply(extract.gt(vcf.SNPs), MARGIN = 1, function(x){ sum(is.na(x)) })
mySNPStats <- data.frame(snp = names(mySNPStats), mn_cov = mySNPStats)

p <- ggplot(mySNPStats) +
  geom_histogram(aes(log(mn_cov+1)))

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_SNPCovHistogram.png"), p, width = 10, height = 10)

# Any populations with all samples with less than 20% missing data
print("These populations have no samples with less than 20% missing data")
unique(mySampleStats$Population)[!unique(mySampleStats$Population)%in%mySampleStats$Population[mySampleStats$per.gt<=20]]

myGT<-extract.gt(vcf.SNPs)
myDP <- extract.gt(vcf.SNPs, element = "DP")
myAD <- as.data.frame(extract.gt(vcf.SNPs, element = "AD"))

myAD.df <- cbind(SNP = rownames(myAD), pivot_longer(myAD, cols = colnames(myAD)))
myAD.df$Ref <- as.numeric(stringr::str_split_i(myAD.df$value, ",", 1))
myAD.df$Alt <- as.numeric(stringr::str_split_i(myAD.df$value, ",", 2))

myAD.df$GT <- pivot_longer(as.data.frame(myGT), cols = colnames(myGT))$value
### 
### ggplot(myAD.df[(myAD.df$Alt+myAD.df$Ref)>=5,]) +
###   geom_histogram(aes(Ref/(Alt+Ref), fill = GT, col = GT), alpha = 0.5, position = "dodge") +
###   coord_cartesian(ylim = c(0,50000))
### 
## Filter out samples with less than 0.8 missing SNP calls
names(mySampleStats)
mySampleStats$sample[mySampleStats$per.gt<=20]
# vcf.SNPs <- vcf.SNPs[samples = mySampleStats$sample[mySampleStats$per.gt<=20]]
## Filter sample file
match(samples_data$ID, colnames(vcf.SNPs@gt)[-1])
samples_data <- samples_data[samples_data$ID%in%(colnames(vcf.SNPs@gt)[-1]), ]

print("Do any samples names not line up with (TRUE is good)")
any(samples_data$ID!=(colnames(vcf.SNPs@gt)[-1]))

print("Check dimensions of sample data and vcf")
dim(samples_data)
dim(vcf.SNPs)

######################################
##### PCA for all samples ####
######################################
vcf.SNPs

## Remove multiallelic snps and snps that are nolonger polymorphic
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs,na.omit = T),]
dim(vcf.SNPs)

## Convert to geno object
print("Converting vcf to geno")
geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)
geno.mat[1:10,1:10]
dim(geno.mat)
table(geno.mat)
geno.mat[geno.mat=="1/1"] <- 2
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="0/0"] <- 0

# Check none of the SNPs are entirely heterozgous and remove them if they are
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}

## Run MDS
dc <- dist(geno.mat)
mds <- cmdscale(dc, k = 4)      

# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9

geno.df <- data.frame(t(geno.mat))
dim(geno.df)

print("Writing out geno file.")
write.table(x = geno.df, file = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read back in geno object
geno <- read.geno(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,".geno"))
dim(geno)

#Conduct PCA
print("Conducting PCA")
geno2lfmm(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,".geno"), 
          paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,".lfmm"), force = TRUE)
#PCA
print("Reading back in PCA")
pc <- pca(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,".lfmm"), scale = TRUE)

pc.sum <- summary(pc)
# Links PCA data to 
print("Merge PCA data with sample data")
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(pc.sum[2,1:6]*100, 1), "%)", sep = "")
pca.comp$sample <- colnames(vcf.SNPs@gt)[-1]
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

print("Creating PCA plots")

pca12.plot <- ggplot(pca.comp) +
  geom_point(aes(pca1, pca2, col = Waterbody)) +
  labs(x = pca.labs[1], y = pca.labs[2])
pca23.plot <- ggplot(pca.comp) +
  geom_point(aes(pca2, pca3, col = Waterbody)) +
  labs(x = pca.labs[2], y= pca.labs[3])
pca45.plot <- ggplot(pca.comp) +
  geom_point(aes(pca4, pca5, col = Waterbody)) +
  labs(x = pca.labs[4], y= pca.labs[5])
pca56.plot <- ggplot(pca.comp) +
  geom_point(aes(pca5, pca6, col = Waterbody)) +
  labs(x = pca.labs[5], y= pca.labs[6])

pca.all.plot <- (pca12.plot + pca23.plot)/(pca45.plot + pca56.plot)

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_PCA.pdf"), pca.all.plot, width = 10, height = 8)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_PCA.png"), pca.all.plot, width = 10, height = 8)

## Add MDS columns
pca.comp$MDS1 <- mds[,1]
pca.comp$MDS2 <- mds[,2]
pca.comp$MDS3 <- mds[,3]

## MDS plots
print("Creating MDS plots")

mds12.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS1, MDS2, col = Waterbody)) +
  labs(x = "MDS1", y = "MDS2")

mds23.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS2, MDS3, col = Waterbody)) +
  labs(x = "MDS2", y = "MDS3")

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS.pdf"), mds12.plot +mds23.plot, width = 10, height = 8)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS.png"), mds12.plot +mds23.plot, width = 10, height = 8)

######################################
##### PCA for paired populations ####
######################################

print("Conducting second PCA just using paired PCA")
# Remove non-species pair locations
paired_sp_waterbodies <- c("DUIN", "LUIB", "CLAC", "OBSE")
# Code to retain only certain samples
vcf.SNPs <- vcf.SNPs[samples = samples_data$ID[samples_data$Waterbody%in%paired_sp_waterbodies]]
## vcf.SNPs
## 
## ## Remove multiallelic snps and snps that are nolonger polymorphic
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs,na.omit = T),]
dim(vcf.SNPs)

## Convert to geno object
print("Converting vcf to geno")
geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)
geno.mat[1:10,1:10]
dim(geno.mat)
table(geno.mat)
geno.mat[geno.mat=="1/1"] <- 2
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="0/0"] <- 0

# Check none of the SNPs are entirely heterozgous and remove them if they are
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}

## Run MDS
dc <- dist(geno.mat)
mds <- cmdscale(dc, k = 6)    

# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9

geno.df <- data.frame(t(geno.mat))
dim(geno.df)

print("Writing out geno file.")
write.table(x = geno.df, file = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read back in geno object
geno <- read.geno(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.geno"))
dim(geno)

#Conduct PCA
geno2lfmm(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.geno"), 
          paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.lfmm"), scale = TRUE)

pc.sum <- summary(pc)
# Links PCA data to 
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(pc.sum[2,1:6]*100, 1), "%)", sep = "")
pca.comp$sample <- colnames(vcf.SNPs@gt)[-1]
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

## Create PCA combined plot
pca12.plot <- ggplot(pca.comp) +
  geom_point(aes(pca1, pca2, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[1], y = pca.labs[2])
pca23.plot <- ggplot(pca.comp) +
  geom_point(aes(pca2, pca3, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[2], y= pca.labs[3])
pca45.plot <- ggplot(pca.comp) +
  geom_point(aes(pca4, pca5, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[4], y= pca.labs[5])
pca56.plot <- ggplot(pca.comp) +
  geom_point(aes(pca5, pca6, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[5], y= pca.labs[6])

pca.all.plot <- (pca12.plot + pca23.plot)/(pca45.plot + pca56.plot) + plot_layout(guides = "collect")

print("Saving PCA plot")
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.pdf"), pca.all.plot, width = 12, height = 8)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.png"), pca.all.plot, width = 12, height = 8)

## PCA strip text
pca.comp.long <- pivot_longer(pca.comp, cols = colnames(pca.comp)[grep("pca", colnames(pca.comp))],
  values_to = "pca_val", names_to = "pca_axis", names_prefix = "pca")

print("Creating PCA strip plot")
pca_strip_plot <- ggplot(pca.comp.long) +
  geom_jitter(aes(pca_val, as.numeric(pca_axis), col = Waterbody, shape = Ecotype), width = 0, height = 0.3) +
  # scale_y_reverse(breaks = 1:6, name =  "PCA axis") +
  theme_bw()

print("Saving pca strip plot")
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA_strip.pdf"), pca_strip_plot)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA_strip.png"), pca_strip_plot)

## Add MDS columns
pca.comp$MDS1 <- mds[,1]
pca.comp$MDS2 <- mds[,2]
pca.comp$MDS3 <- mds[,3]
pca.comp$MDS4 <- mds[,4]
pca.comp$MDS5 <- mds[,5]
pca.comp$MDS6 <- mds[,6]

## MDS plots
print("Creating MDS plots")

mds12.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS1, MDS2, col = Waterbody, shape = Ecotype)) +
  geom_text_repel(data = pca.comp[pca.comp$sample=="Uist22CLAM4",], aes(MDS1, MDS2, label = sample),
                   alpha = 0.8, nudge_y = 10, min.segment.length = 0)
  labs(x = "MDS1", y = "MDS2")

mds23.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS1, MDS3, col = Waterbody, shape = Ecotype)) +
  geom_text_repel(data = pca.comp[pca.comp$sample=="Uist22CLAM4",], aes(MDS1, MDS3, label = sample),
                   alpha = 0.8, nudge_x = -20, min.segment.length = 0)
  labs(x = "MDS2", y = "MDS3")

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS_paired.pdf"), mds12.plot + mds23.plot, width = 10, height = 6)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS_paired.png"), mds12.plot + mds23.plot, width = 10, height = 6)

## PCA strip text
mds.comp.long <- pivot_longer(pca.comp, cols = colnames(pca.comp)[grep("MDS", colnames(pca.comp))],
values_to = "mds_val", names_to = "mds_axis", names_prefix = "MDS")

mds_strip_plot <- ggplot(mds.comp.long) +
  geom_jitter(aes(mds_val, as.numeric(mds_axis), col = Waterbody, shape = Ecotype), width = 0, height = 0.3) +
  scale_y_reverse(breaks = 1:6, name =  "MDS axis") +
  theme_bw()

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_MDS_strip.pdf"), mds_strip_plot)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_MDS_strip.png"), mds_strip_plot)

# # # # # # # # # # # # # # # #
####### Allele overlap ######
# # # # # # # # # # # # # # # #
# install.packages("ggVennDiagram")
library(ggVennDiagram)

#Read back in geno object
geno <- read.geno(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.geno"))
geno[geno==9] <- NA

## Make sure geno is polarized (minor allele is equal to 2)
allele.polz <- function(x){
  minor.allele.freq <- sum(x[(x==1|x==2)&!is.na(x)]/2)/sum(as.numeric(!is.na(x))) # Calculate 2 allele frequency
  if(minor.allele.freq > 0.5){ # If allele freq is greater than 0.5
    re.polz <- as.integer(factor(x, levels = c(2,1,0)))-1 ## Invert allele symbols (2 to 0, 0 to 2)
  } else { re.polz <- x } # Else fo nothing
  return(re.polz)
}

# Make allele "2" the minor allele
geno.polz <- apply(geno, MARGIN = 2, allele.polz)

# Check allele frequency
minor.allele.freq <- apply(geno.polz, MARGIN = 2, function(x) {
  sum(x[(x==1|x==2)&!is.na(x)]/2)/sum(as.numeric(!is.na(x)))
})

# Get populations
pops <- pca.comp$Population
nsnps <- ncol(geno) # Get number of snps
nind <- nrow(geno) # Get number of samples
# Create unique ID for each snp
snps.names <- paste0("snp", 1:nsnps)

# Create list which contains a string of all the snps foud within each population
## All Anad populations
print("Starting Private allele calcs")
venn.data.anad <- list()
for(pop in c("CLAM", "DUIM", "OBSM", "LUIM")){
  # Extract name of each snp that is foud in each population
  venn.data.anad[[pop]] <- snps.names[apply(geno.polz[pops==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}
# Create Venn diagram
p.venn.anad <- ggVennDiagram(venn.data.anad) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))
print("1")
# All resident popuation
venn.data.resi <- list()
for(pop in c("CLAC", "DUIN", "OBSE", "LUIB")){
  venn.data.resi[[pop]] <- snps.names[apply(geno.polz[pops==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}
## Plot Venn diagram
p.venn.resi <- ggVennDiagram(venn.data.resi) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))
print("2")
# Resi vs Eco popuation
venn.data.eco <- list()
Ecos <- pca.comp$Ecotype
for(pop in c(c("anad", "resi"))){
  venn.data.eco[[pop]] <- snps.names[apply(geno.polz[Ecos==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}

## Plot Venn diagram
p.venn.eco <- ggVennDiagram(venn.data.eco) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))
print("3")
# All resident popuation and all combined anad
venn.data.resi.p.tand <- list()
pops <- pca.comp$Population
pops[pops%in%c("CLAM", "DUIM", "OBSM", "LUIM")] <- "anad"
for(pop in unique(pops)){
  venn.data.resi.p.tand[[pop]] <- snps.names[apply(geno.polz[pops==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}
## Plot Venn diagram
p.venn.resi.p.tand <- ggVennDiagram(venn.data.resi.p.tand) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))
print("4")
# All resident popuation and all combined anad
venn.data.anad.p.tresi <- list()
pops <- pca.comp$Population
pops[pops%in%c("CLAC", "DUIN", "OBSE", "LUIB")] <- "resi"
for(pop in unique(pops)){
  venn.data.anad.p.tresi[[pop]] <- snps.names[apply(geno.polz[pops==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}
## Plot Venn diagram
p.venn.anad.p.tresi <- ggVennDiagram(venn.data.anad.p.tresi) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))

print("5")
# All resident popuation and all combined anad
venn.data.waterbody<- list()
pops <- pca.comp$Waterbody
for(pop in unique(pops)){
  venn.data.waterbody[[pop]] <- snps.names[apply(geno.polz[pops==pop,], MARGIN = 2, function(x) any(x==1|x==2))]
}
## Plot Venn diagram
p.venn.waterbody <- ggVennDiagram(venn.data.waterbody) + scale_fill_gradient(low="grey90",high = "red") + theme(plot.background = element_rect(fill = "white"))


# Save
print("Saving private allele figures")
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_resi.png"), p.venn.resi)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_anad.png"), p.venn.anad)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_Ecotype.png"), p.venn.eco)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_resi-pops_all-anad.png"), p.venn.resi.p.tand )
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_anad-pops_all-resi.png"), p.venn.anad.p.tresi )
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MA_private_alleles_waterbody.png"), p.venn.waterbody )


###### Rarefaction of private alleles #######

# Convert vcf to genind
geno.genind <- vcfR2genind(vcf.SNPs, return.alleles = F)
# Assign populations
rownames(geno.genind@tab)==pca.comp$sample
geno.genind@pop <- as.factor(pca.comp$Population)

## All samples number of private alleles
pa_raw <- private_alleles(geno.genind, count.alleles = F)
pa_raw_total <- rowMeans(pa_raw)

## What min number of samples from a population
min_n <- min(table(geno.genind@pop))

# Create 100 permuations of samples selection
sample_perms <- lapply(1:100, function(x) 
  samples <- sort(unlist(lapply(unique(geno.genind@pop), function(p) {
    sample(which(geno.genind@pop==p), min_n)
  })))
)
## Remove any duplicated sample selection
sample_perms <- sample_perms[!duplicated(unlist(lapply(sample_perms, function(x) paste0(x, collapse = ""))))]

## Calculate number of private alleles when subsetting to lower sample sets of each population
priv_table <- lapply(sample_perms, function(x) {
  priv_df <- private_alleles(geno.genind[x,], count.alleles = F)
  return(rowSums(priv_df))
  }
)
## Convert to tidy format
priv_dataframe <- pivot_longer(data.frame(do.call("rbind", priv_table)), cols = unique(geno.genind@pop), names_to = "Population", values_to = "Private_Alleles")

# Reorder Ppulations to be grouped by Ecotype
priv_dataframe$Population <- factor(priv_dataframe$Population, levels = c("CLAC", "DUIN", "LUIB", "OBSE", 
                                                                          "CLAM", "DUIM", "LUIM", "OBSM"))

# Convert to percentage
priv_dataframe$Private_Alleles_pcent <- (priv_dataframe$Private_Alleles/dim(geno.genind@tab)[2])*100

# Summarise rarefaction
pa.summary <- priv_dataframe %>%
group_by(Population) %>%
  summarise(mn = mean(Private_Alleles_pcent), sd.pa = sd(Private_Alleles_pcent))
# Plot distibution of private alleles
priv_plot <- ggplot(priv_dataframe) +
  geom_jitter(aes(Population, Private_Alleles_pcent), width = 0.2, height = 0, col = "grey50") +
  geom_point(data = pa.summary, aes(Population, mn), col = "black", size = 3) +
  geom_segment(data = pa.summary, aes(Population, y = mn+sd.pa, yend = mn-sd.pa), col = "black", linewidth = 2, lineend = "round") +
  theme_bw() +
  lims(y = c(0, max(priv_dataframe$Private_Alleles_pcent))) +
  labs(x = "Population", y = "Number of Private Alleles (%)")

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_Rarefaction_private_alleles_n", min_n,".png"), priv_plot)

### nj tree plot ###
### ## Not currently working on Ada as ggtree is not installing
### 
### #  Calc nj 
### nj.data <- nj(dist(tab(geno.genind, freq=TRUE)))
### 
### # Create tree plot
### plot.tree <- ggtree(nj.data, layout = "daylight")
### #  Combine with sample data
### plot.tree <- plot.tree %<+% pca.comp
### 
### ## Custom tip colours
### plot.tree <- plot.tree + geom_tippoint(aes(color=Ecotype, fill  = Waterbody), size=3, shape  = 21, stroke = 2) +
###   scale_color_manual(values = c("black", "grey"))
### 
### # Save
### ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_njtree.png"), plot.tree)
### 

# # # # # # # # # # # # # # # #
####### Kinship analysis ######
# # # # # # # # # # # # # # # #

#install.packages("popkin")
library(popkin)
dir.create(paste0(plot.dir, "/kinship/"))

# Load geno file and change missing values to "NA
geno.kin <- read.geno(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.geno"))
geno.kin[geno.kin=="9"] <- NA

# Choose subpopulations
#Check dimensions are correct
paired_samples <- samples_data[match(colnames(vcf.SNPs@gt)[-1], samples_data$ID),]

dim(paired_samples)[1]==dim(geno.kin)[1]
subpops <- paired_samples$Population
# Reorder samples
kin.plot.order <- order(paste(paired_samples$Ecotype ,paired_samples$Population))
paired_samples <- paired_samples[kin.plot.order,]
geno.kin <- geno.kin[kin.plot.order,]
subpops <- subpops[kin.plot.order]
subpops.site.sub <- paired_samples$Population[kin.plot.order]

# plot as png
kinship <- popkin(t(geno.kin), subpops = subpops)
png(paste0(plot.dir, "/kinship/Kinship_popkin_baseplot_", SNP.library.name,".png"), width = 2000, height = 2000)
plot_popkin(
  kinship,
  labs = subpops,
  # shared bottom and left margin value, to make space for labels
  mar = 1
)
dev.off()

## Plot as pdf
pdf(paste0(plot.dir, "/kinship/Kinship_popkin_baseplot_", SNP.library.name,".pdf"), width = 20, height = 20)
plot_popkin(
  kinship,
  labs = subpops,
  # shared bottom and left margin value, to make space for labels
  mar = 1
)
dev.off()

## Convert to tidy format
# rename rows
colnames(kinship) <- paired_samples$ID
rownames(kinship) <- paired_samples$ID
# convert to dataframe
kinship.df <- as.data.frame(kinship)
# Add first sample compares
kinship.df$ID.1 <- rownames(kinship.df)
# Pivot to longer
kinship.df <- pivot_longer(kinship.df, cols = paired_samples$ID, values_to = "Kinship", names_to = "ID.2")

# Remove data comparing self with self
# remove Sample Obsm_491
any(!(kinship.df$ID.1=="Obsm_641"|kinship.df$ID.2=="Obsm_641"))
kinship.df <-  kinship.df[!(kinship.df$ID.1=="Obsm_641"|kinship.df$ID.2=="Obsm_641"),]

## ad pops 
# Make order to take into account Ecotype
kinship.df$Population_1 <- factor(paired_samples$Population[match(kinship.df$ID.1, paired_samples$ID)], levels = c("CLAC", "DUIN", "LUIB", "OBSE", 
                                                                                                                          "CLAM", "DUIM", "LUIM", "OBSM"))
kinship.df$Population_2 <- factor(paired_samples$Population[match(kinship.df$ID.2, paired_samples$ID)], levels = c("CLAC", "DUIN", "LUIB", "OBSE", 
                                                                                                                           "CLAM", "DUIM", "LUIM", "OBSM"))
paired_samples$Population <- factor(paired_samples$Population, levels = c("CLAC", "DUIN", "LUIB", "OBSE", 
                                             "CLAM", "DUIM", "LUIM", "OBSM"))

# Set all kinship calc comparing to self equal to NA
kinship.df$Kinship[kinship.df$ID.1==kinship.df$ID.2] <- NA
## Order samples to be same as population
kinship.df$ID.1 <- factor(kinship.df$ID.1, levels = paired_samples$ID[order(paired_samples$Population)])
kinship.df$ID.2 <- factor(kinship.df$ID.2, levels = paired_samples$ID[order(paired_samples$Population)])

## Create annotation data
## df of counts
pop_bar <- kinship.df %>%
  mutate(x = ID.1, y = 1)

# Get population sizes in the order of samples
pop_sizes <- paired_samples %>%
  arrange(Population) %>%
  count(Population)

# Compute cumulative positions for line placement
breaks <- cumsum(pop_sizes$n) + 0.5

# Calculate the midpoint of each population block
pop_sizes <- pop_sizes %>%
  mutate(start = cumsum(lag(n, default = 0)) + 1,
         end = cumsum(n),
         center = (start + end) / 2)

## Plot
kin.plot <- ggplot(kinship.df) +
  geom_tile(aes(ID.1, ID.2, fill = Kinship)) +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_vline(xintercept = breaks) +
  geom_hline(yintercept = breaks) +
  annotate("text", x = pop_sizes$center, y = 1,
           label = pop_sizes$Population, angle = 0, vjust = 2, size = 3.5) +
  
  annotate("text", y = pop_sizes$center, x = 1,
           label = pop_sizes$Population, angle = 90, vjust = -1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(20, 20, 20, 30))

ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_", SNP.library.name,".png"), kin.plot, width = 7, height = 6)
ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_", SNP.library.name,".pdf"), kin.plot, width = 7, height = 6)

############################
 ##### Nj dist plot #####
############################



#Calculates structure for samples from K=1 to k=15
max.K <- 6
# MAY NEED TO PAUSE ONEDRIVE
# File names are becoming too Long
obj.at <- snmf(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.geno"), K = 1:max.K, ploidy = 2, entropy = T,
             CPU = 12, project = "new", repetitions = 10, alpha = 100)
stickleback.snmf <- load.snmfProject(file = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.snmfProject"))
stickleback.snmf.sum <- summary(stickleback.snmf)

plot(stickleback.snmf, col = "blue4", cex = 1.4, pch = 19)

ce <- cbind(1:max.K, t(stickleback.snmf.sum$crossEntropy))
colnames(ce) <- c("K", "min","mean","max")
ce <- data.frame(ce)

summary(stickleback.snmf)

k.min <- which.min(ce$mean)
#Choose K
if(k.min==1) {
  print("Min K equal to 1, setting K to 2")
  K <- 2
  } else {
    print(paste("Min K was greater than 1, setting K to", k.min))
    K <- k.min
    }

# K <- 4
ce.plot <- ggplot(ce) +
  geom_point(aes(K, mean), size = 2, shape = 19) +
  geom_errorbar(aes(x = K, ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +
  ggtitle(paste("Minimum mean cross-entropy value K =", which.min(ce$mean))) +
  ylab("Cross-entropy") +
  theme_bw()
ce.plot
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_LEA_K1-",max.K,"_cross_entropy.pdf"), plot = ce.plot)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_LEA_K1-",max.K,"_cross_entropy.jpg"), plot = ce.plot)

best <- which.min(cross.entropy(stickleback.snmf, K = K))

qmatrix = Q(stickleback.snmf, K = K, run = best)
dim(qmatrix)
# Tidy data for plotting
qtable <-  cbind(colnames(vcf.SNPs@gt)[-1], rep(1:K, each = dim(qmatrix)[1]), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","Qid", "Q")

qtable$Q <- as.numeric(qtable$Q)
qtable <- merge(qtable, samples_data, by.x = "sample", by.y = "ID")
#qtable$sample_Lat <- paste(qtable$Lat, qtable$sample, sep = "_")

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")
cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")

head(qtable)
v <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = as.factor(Qid)), position = "stack", width = 1, col = "black") +
  scale_fill_manual(values = cbPalette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  facet_grid(~Ecotype+Population, drop = T, scales = "free", space = "free") +
  ylab(label = paste("K =", K))
v

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_LEA_K",K,".pdf"), v, width = 18, height = 6)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_LEA_K",K,".png"), v, width = 18, height = 6)


#Creates multiple K barcharts

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")

s <- list()

max.K <- min(c(max.K, 6))

for(i in 2:max.K){
  best <- which.min(cross.entropy(stickleback.snmf, K = i))
  qmatrix = Q(stickleback.snmf, K = i, run = best)
  dim(qmatrix)
  # Tidy data for plotting
  qtable <-  cbind(colnames(vcf.SNPs@gt)[-1], rep(1:i, each = dim(qmatrix)[1]), c(qmatrix[,1:i]))
  qtable <-  data.frame(qtable)
  colnames(qtable) <- c("sample","Qid", "Q")
  
  qtable$Q <- as.numeric(qtable$Q)
  qtable <- merge(qtable, samples_data, by.x = "sample", by.y = "ID")
    
  s[[i]] <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = as.factor(Qid)), position = "stack", width = 1, col = "black") +
  scale_fill_manual(values = cbPalette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  facet_grid(~Ecotype+Population, drop = T, scales = "free", space = "free") +
  ylab(label = paste("K =", i))
  
}

plot2 <- s[[2]]
for(i in 3:max.K){
  plot2 <- plot2 / s[[i]]
}
s <- plot2

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_LEA_barplot_1-",max.K,".pdf"), plot=s, height=20, width=15)
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_LEA_barplot_1-",max.K,".jpg"), plot=s, height=20, width=15)

