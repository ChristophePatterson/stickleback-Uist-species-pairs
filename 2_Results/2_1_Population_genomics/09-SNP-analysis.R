## Ada LEA and PCR plots
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
library(ggnewscale)

#Not currently installed
#library(treedataverse)

## Colorblind palette
cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
geno.file <- args[1]
# geno.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.geno"

vcf.ver <- args[2]
# vcf.ver <- "GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20"
# Remove file extension
SNP.library.name <- basename(gsub(".geno", "", geno.file))

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
geno <- read.geno(geno.file)
dim(geno)
# Get an read sample information
samples_data <- data.frame(ID = readLines(paste0(geno.file, ".samples")))
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)

colnames(geno) <- samples_data$ID
## Fill in missing data
samples_data$Ecotype[is.na(samples_data$Ecotype)] <- "Unknown"

### Missing genotype assesment
# Stats
mySampleStats <- apply(geno, MARGIN = 2, function(x){ sum(x == 9) })
mySampleStats <- data.frame(sample = names(mySampleStats), per.gt = mySampleStats/nrow(geno)*100)
dim(mySampleStats)
mySampleStats <- merge(mySampleStats, samples_data[,c(1,7,8,9,10,11,12,13,14)], by.x = "sample", by.y = "ID")

summary(mySampleStats)
p <- ggplot(mySampleStats) +
  geom_boxplot(aes(x = Population, y = per.gt, col = Population), outlier.colour = NA) +
  geom_jitter(aes(x = Population, y = per.gt, col = Population), height = 0, width = 0.2) 

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_missGT.png"), p, width = 10, height = 10)

mySNPStats <- data.frame(SNPmiss = apply(geno, MARGIN = 1, function(x){ sum(x == 9 ) }))
mySNPStats <- data.frame(snp = names(mySNPStats), mn_cov = mySNPStats)

p <- ggplot(mySNPStats) +
  geom_histogram(aes(SNPmiss))

ggsave("test.png", p, width = 10, height = 10)

ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_SNPCovHistogram.png"), p, width = 10, height = 10)

# Any populations with all samples with less than 20% missing data
print("These populations have no samples with less than 40% missing data")
unique(mySampleStats$Population)[!unique(mySampleStats$Population)%in%mySampleStats$Population[mySampleStats$per.gt<=40]]

dim(geno)
# Filter out samples with greater than 40 missing GT
geno <- geno[,mySampleStats$sample[mySampleStats$per.gt<=40]]
## Filter sample file
samples_data <- samples_data[samples_data$ID%in%(mySampleStats$sample[mySampleStats$per.gt<=40]),]

no.longer.poly <- apply(geno, MARGIN = 1, function(x) length(unique(x[x!=9]))>1)
geno <- geno[no.longer.poly,]

print("Do any samples names not line up with (FALSE is good)")
any(samples_data$ID!=(colnames(geno)))

## remove any SNPs that are nolonger polymorphic

print("Check dimensions of sample data and geno")
dim(samples_data)
dim(geno)

######################################
##### PCA for all samples ####
######################################

write.geno(t(geno), paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,".geno"))

## Run MDS
geno.mat <- t(geno)
geno.mat[geno.mat==9] <- NA
dc <- dist(geno.mat)
mds <- cmdscale(dc, k = 4)      

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
pca.comp$sample <- colnames(geno)
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

print("Creating PCA plots")

pca12.plot <- ggplot(pca.comp) +
  geom_point(aes(pca1, pca2, col = Waterbody)) +
  labs(x = pca.labs[1], y = pca.labs[2]) + theme_bw()
pca23.plot <- ggplot(pca.comp) +
  geom_point(aes(pca2, pca3, col = Waterbody)) +
  labs(x = pca.labs[2], y= pca.labs[3]) + theme_bw()
pca45.plot <- ggplot(pca.comp) +
  geom_point(aes(pca4, pca5, col = Waterbody)) +
  labs(x = pca.labs[4], y= pca.labs[5]) + theme_bw()
pca56.plot <- ggplot(pca.comp) +
  geom_point(aes(pca5, pca6, col = Waterbody)) +
  labs(x = pca.labs[5], y= pca.labs[6]) + theme_bw()

pca.all.plot <- (pca12.plot + pca23.plot)/(pca45.plot + pca56.plot) + plot_layout(guides = 'collect')

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
  labs(x = "MDS1", y = "MDS2") + theme_bw()

mds23.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS2, MDS3, col = Waterbody)) +
  labs(x = "MDS2", y = "MDS3") + theme_bw()

mdsplot <- (mds12.plot + mds23.plot) + plot_layout(guides = 'collect')

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS.pdf"), mdsplot, width = 10, height = 6)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS.png"), mdsplot, width = 10, height = 6)

######################################
##### PCA for paired populations ####
######################################

print("Conducting second PCA just using paired PCA")
# Remove non-species pair locations
paired_sp_waterbodies <- c("DUIN", "LUIB", "CLAC", "OBSE")
# Code to retain only certain samples
geno <- geno[,samples_data$ID[samples_data$Waterbody%in%paired_sp_waterbodies]]
## Filter sample file
samples_data <- samples_data[samples_data$Waterbody%in%paired_sp_waterbodies,]
any(colnames(geno)!=samples_data$ID)


## ## Remove multiallelic snps and snps that are nolonger polymorphic
no.longer.poly <- apply(geno, MARGIN = 1, function(x) length(unique(x[x!=9]))>1)
geno <- geno[no.longer.poly,]

#Read back in geno object
write.geno(t(geno), paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.geno"))

## Run MDS
geno.mat <- t(geno)
geno.mat[geno.mat==9] <- NA
dc <- dist(geno.mat)
mds <- cmdscale(dc, k = 6)   

#Conduct PCA
geno2lfmm(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.geno"), 
          paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.lfmm"), scale = TRUE)
pc.sum <- summary(pc)

# Links PCA data to 
print("Merge PCA data with sample data")
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(pc.sum[2,1:6]*100, 1), "%)", sep = "")
pca.comp$sample <- colnames(geno)
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

## Create PCA combined plot
pca12.plot <- ggplot(pca.comp) +
  geom_point(aes(pca1, pca2, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[1], y = pca.labs[2]) +
  scale_color_manual(values = cbPalette) + theme_bw()
pca23.plot <- ggplot(pca.comp) +
  geom_point(aes(pca2, pca3, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[2], y= pca.labs[3]) +
  scale_color_manual(values = cbPalette) + theme_bw()
pca45.plot <- ggplot(pca.comp) +
  geom_point(aes(pca4, pca5, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[4], y= pca.labs[5]) +
  scale_color_manual(values = cbPalette) + theme_bw()
pca56.plot <- ggplot(pca.comp) +
  geom_point(aes(pca5, pca6, col = Waterbody, shape = Ecotype), size = 3) +
  labs(x = pca.labs[5], y= pca.labs[6]) +
  scale_color_manual(values = cbPalette) + theme_bw()

pca.all.plot <- (pca12.plot + pca23.plot)/(pca45.plot + pca56.plot) + plot_layout(guides = "collect")
pca.12.plot <- (pca12.plot + pca23.plot) + plot_layout(guides = "collect")

print("Saving PCA plot")
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.pdf"), pca.all.plot, width = 10, height = 8)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.png"), pca.all.plot, width = 10, height = 8)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.pdf"), pca.12.plot, width = 10, height = 5)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired_PCA.png"), pca.12.plot, width = 10, height = 5)

## PCA strip text
pca.comp.long <- pivot_longer(pca.comp, cols = colnames(pca.comp)[grep("pca", colnames(pca.comp))],
  values_to = "pca_val", names_to = "pca_axis", names_prefix = "pca")

print("Creating PCA strip plot")
pca_strip_plot <- ggplot(pca.comp.long) +
  geom_jitter(aes(pca_val, as.numeric(pca_axis), col = Waterbody, shape = Ecotype), width = 0, height = 0.3) +
  scale_y_reverse(breaks = 1:6, name =  "PCA axis") +
  labs(x = "PCA axis value") +
  scale_color_manual(values = cbPalette) + theme_bw()

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
                   alpha = 0.8, nudge_y = 10, min.segment.length = 0) +
  scale_color_manual(values = cbPalette) +
  labs(x = "MDS1", y = "MDS2") + theme_bw()

mds23.plot <- ggplot(pca.comp) +
  geom_point(aes(MDS1, MDS3, col = Waterbody, shape = Ecotype)) +
  geom_text_repel(data = pca.comp[pca.comp$sample=="Uist22CLAM4",], aes(MDS1, MDS3, label = sample),
                   alpha = 0.8, nudge_x = -20, min.segment.length = 0) +
  scale_color_manual(values = cbPalette) +
  labs(x = "MDS2", y = "MDS3") + theme_bw()

mdsplot <- (mds12.plot + mds23.plot) + plot_layout(guides = 'collect')

ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS_paired.pdf"), width = 8, height = 4)
ggsave(filename = paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_MDS_paired.png"), width = 8, height = 4)

## PCA strip text
mds.comp.long <- pivot_longer(pca.comp, cols = colnames(pca.comp)[grep("MDS", colnames(pca.comp))],
values_to = "mds_val", names_to = "mds_axis", names_prefix = "MDS")

mds_strip_plot <- ggplot(mds.comp.long) +
  geom_jitter(aes(mds_val, as.numeric(mds_axis), col = Waterbody, shape = Ecotype), width = 0, height = 0.3) +
  scale_y_reverse(breaks = 1:6, name =  "MDS axis") +
  scale_color_manual(values = cbPalette) +
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
print("starting private alleles calculation")
# Convert geno to genind
dim(geno)
geno.genind <- as.genind(geno)

# Assign populations
rownames(geno.genind@tab) <- pca.comp$sample
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

subset_snps <- round(ncol(geno.genind@tab)*0.1)

## Calculate number of private alleles when subsetting to lower sample sets of each population
print("Starting private allele permutations")
priv_table <- lapply(sample_perms, function(x) {
  priv_df <- private_alleles(geno.genind[x,sample(1:ncol(geno.genind@tab),subset_snps)], count.alleles = F)
  return(rowSums(priv_df))
  }
)
## Convert to tidy format
priv_dataframe <- pivot_longer(data.frame(do.call("rbind", priv_table)), cols = unique(geno.genind@pop), names_to = "Population", values_to = "Private_Alleles")

# Reorder Ppulations to be grouped by Ecotype
priv_dataframe$Population <- factor(priv_dataframe$Population, levels = c("CLAC", "DUIN", "LUIB", "OBSE", 
                                                                          "CLAM", "DUIM", "LUIM", "OBSM"))

# Convert to percentage
priv_dataframe$Private_Alleles_pcent <- (priv_dataframe$Private_Alleles/subset_snps)*100

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


SNP.library.name
ggsave(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name, "_Rarefaction_private_alleles_n", min_n[1],".png"), priv_plot)

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
print("starting kinship analysis")
dir.create(paste0(plot.dir, "/kinship/"))

# Load geno file and change missing values to "NA
geno.kin <- read.geno(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/",SNP.library.name,"_paired.geno"))
geno.kin[geno.kin=="9"] <- NA

# Choose subpopulations
#Check dimensions are correct
paired_samples <- samples_data[samples_data$Waterbody%in%paired_sp_waterbodies,]

dim(paired_samples)[1]==dim(geno.kin)[1]
subpops <- paired_samples$Population
# Reorder samples
kin.plot.order <- order(paste(paired_samples$Ecotype ,paired_samples$Population))
paired_samples <- paired_samples[kin.plot.order,]
geno.kin <- geno.kin[kin.plot.order,]
subpops <- subpops[kin.plot.order]
subpops.site.sub <- paired_samples$Population[kin.plot.order]

# Calculate popkin and replace diagonals with inbreeding
kinship <- inbr_diag(popkin(t(geno.kin), subpops = subpops))

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

## Order samples to be same as population
kinship.df$ID.1 <- factor(kinship.df$ID.1, levels = paired_samples$ID[order(paired_samples$Population)])
kinship.df$ID.2 <- factor(kinship.df$ID.2, levels = paired_samples$ID[order(paired_samples$Population)])

# Set all kinship calc comparing to self equal to NA
inbr.df <- kinship.df[kinship.df$ID.1==kinship.df$ID.2,]
kinship.df$Kinship[kinship.df$ID.1==kinship.df$ID.2] <- NA

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
  new_scale_fill() +
  geom_tile(data = inbr.df, aes(ID.1, ID.2, fill = Kinship)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkgreen") +
  geom_vline(xintercept = breaks) +
  geom_hline(yintercept = breaks) +
  annotate("text", x = pop_sizes$center, y = 1,
           label = pop_sizes$Population, angle = 0, vjust = 2, size = 3.5) +
  
  annotate("text", y = pop_sizes$center, x = 1,
           label = pop_sizes$Population, angle = 90, vjust = -1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "Inbreeding") +
  coord_cartesian(clip = "off") +
  theme(text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(20, 20, 20, 30))

ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_", SNP.library.name,".png"), kin.plot, width = 7, height = 6)
ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_", SNP.library.name,".pdf"), kin.plot, width = 7, height = 6)

pops_kin_plot <- ggplot(kinship.df[kinship.df$Population_1==kinship.df$Population_2,], aes(Population_1, Kinship, fill = Population_1)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(height = 0, width = 0.2) +
  scale_fill_manual(values = cbPalette) +
  theme_bw() +
  labs(fill = "Population", x = "Population")


ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_intraPop_", SNP.library.name,".png"), pops_kin_plot, width = 7, height = 6)
ggsave(paste0(plot.dir, "/kinship/Kinship_popkin_ggplot_intraPop_", SNP.library.name,".pdf"), pops_kin_plot, width = 7, height = 6)

pops_inbr_plot <- ggplot(inbr.df, aes(Population_1, Kinship, fill = Population_1)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(height = 0, width = 0.2) +
  scale_fill_manual(values = cbPalette) +
  theme_bw() +
  labs(fill = "Population", x = "Population")

ggsave("test.png", pops_inbr_plot )
ggsave(paste0(plot.dir, "/kinship/Interbreeding_popkin_ggplot_intraPop_", SNP.library.name,".png"), pops_inbr_plot, width = 7, height = 6)
ggsave(paste0(plot.dir, "/kinship/Interbreeding_popkin_ggplot_intraPop_", SNP.library.name,".pdf"), pops_inbr_plot, width = 7, height = 6)

## Most inbreed individuals
geno.genind.inbr <- geno.genind[which.min(inbr.df$Kinship),]

inbr.df[order(inbr.df$Kinship),]
totalSamp <- dim(geno.genind@tab)[1]
totalAll <- dim(geno.genind@tab)[2]
het.df <- data.frame(sample = row.names(geno.genind@tab), hetO = apply(geno.genind@tab, MARGIN = 1, function(x) sum(x[!is.na(x)] == 1)/(totalAll-sum(is.na(x)))),
                      missGT = apply(geno.genind@tab, MARGIN = 1, function(x) sum(is.na(x))/totalAll))


all.pair.stats <- merge(het.df, pca.comp, by = "sample")
all.pair.stats <- merge(all.pair.stats, inbr.df[,c("ID.1", "Kinship")], by.x = "sample", by.y = "ID.1")

plot.pair.stats <- ggplot(all.pair.stats, aes(hetO, MDS1, col = missGT, shape = Waterbody)) +
  geom_point() +
  # scale_fill_manual(values = cbPalette) +
  theme_bw() +
  labs(fill = "Population", x = "het0")

ggsave(paste0("test.png"), plot.pair.stats , width = 7, height = 6)

kinship.df[kinship.df$ID.1=="Uist22633"|kinship.df$ID.2=="Uist22633",][which.max(kinship.df[kinship.df$ID.1=="Uist22633"|kinship.df$ID.2=="Uist22633",]$Kinship),]

all.pair.stats[all.pair.stats$sample=="Uist22609",]

############################
 ##### Nj dist plot #####
############################



#Calculates structure for samples from K=1 to k=15
max.K <- 6
# MAY NEED TO PAUSE ONEDRIVE
# File names are becoming too Long
obj.at <- snmf(paste0(plot.dir, "/LEA_PCA/", SNP.library.name, "/", SNP.library.name,"_paired.geno"), K = 1:max.K, ploidy = 2, entropy = T,
             CPU = 12, project = "new", repetitions = 100, alpha = 100)
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
qtable <-  cbind(row.names(geno.genind@tab), rep(1:K, each = dim(qmatrix)[1]), c(qmatrix[,1:K]))
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
  qtable <-  cbind(row.names(geno.genind@tab), rep(1:i, each = dim(qmatrix)[1]), c(qmatrix[,1:i]))
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

