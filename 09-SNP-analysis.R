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

#Not currently installed
#library(poppr)
#library(ggnewscale)

SNP.library.name <- "stickleback"

# Test whether working on HPC or laptop and set working directory accordingly
# Laptop test
if(grepl(getwd(), pattern = "C:/Users/mbzcp2/")){
  dir.path <-"C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/"
  plot.dir <- "C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/"
}
# HPC test
if(grepl(getwd(), pattern = "/gpfs01/home/mbzcp2")){
  dir.path <-"/gpfs01/home/mbzcp2/data/sticklebacks/"
  plot.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/"
}
## Create directory is not already
dir.create(plot.dir)
dir.create(paste0(plot.dir, "/LEA_PCA/"))

vcf.SNPs <- read.vcfR(paste0(dir.path, "vcfs/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz"),
                      verbose = T)
# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 

colnames(vcf.SNPs@gt)
## Remove one of the dupicalted samples
vcf.SNPs <- vcf.SNPs[samples = colnames(vcf.SNPs@gt)[colnames(vcf.SNPs@gt)!="Obsm_641"]]

# Get an read sample information
samples_data <- data.frame(ID = colnames(vcf.SNPs@gt)[-1])
samples <- read.csv("bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)
samples_data <- samples_data[samples_data$ID%in%colnames(vcf.SNPs@gt),]
samples_data <- samples_data[match(samples_data$ID, (colnames(vcf.SNPs@gt)[-1])),]
print("Do any samples names not line up with (False is good)")
any(!samples_data$ID==(colnames(vcf.SNPs@gt)[-1]))
## Fill in missing data
samples_data$Ecotype[is.na(samples_data$Ecotype)] <- "Unknown"

### Missing genotype assesment
# Stats
myMiss <- apply(extract.gt(vcf.SNPs), MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- data.frame(sample = names(myMiss), per.gt = myMiss/nrow(vcf.SNPs)*100)
dim(myMiss)
myMiss <- merge(myMiss, samples_data[,c(1,7,8,9,10,11,12,13,14)], by.x = "sample", by.y = "ID")

summary(myMiss)
p <- ggplot(myMiss) +
  geom_boxplot(aes(x = Population, y = per.gt, col = Population), outlier.colour = NA) +
  geom_jitter(aes(x = Population, y = per.gt, col = Population), height = 0, width = 0.2) 

ggsave("test.png", p, width = 10, height = 10)
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
# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9

geno.df <- data.frame(t(geno.mat))
dim(geno.df)

print("Writing out geno file.")
write.table(x = geno.df, file = paste0(plot.dir,SNP.library.name,".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read back in geno object
geno <- read.geno(paste0(plot.dir,SNP.library.name,".geno"))
dim(geno)

#Conduct PCA
geno2lfmm(paste0(plot.dir, SNP.library.name,".geno"), 
          paste0(plot.dir,SNP.library.name,".lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(plot.dir,SNP.library.name,".lfmm"), scale = TRUE)

pc.sum <- summary(pc)
# Links PCA data to 
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(pc.sum[2,1:6]*100, 1), "%)", sep = "")
pca.comp$sample <- colnames(vcf.SNPs@gt)[-1]
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

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

ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA.pdf"), pca.all.plot, width = 15, height = 15)
ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA.png"), pca.all.plot, width = 15, height = 15)

######################################
##### PCA for paired populations ####
######################################

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
# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9

geno.df <- data.frame(t(geno.mat))
dim(geno.df)

print("Writing out geno file.")
write.table(x = geno.df, file = paste0(plot.dir,SNP.library.name,"_paired.geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read back in geno object
geno <- read.geno(paste0(plot.dir,SNP.library.name,"_paired.geno"))
dim(geno)

#Conduct PCA
geno2lfmm(paste0(plot.dir, SNP.library.name,"_paired.geno"), 
          paste0(plot.dir,SNP.library.name,"_paired.lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(plot.dir,SNP.library.name,"_paired.lfmm"), scale = TRUE)

pc.sum <- summary(pc)
# Links PCA data to 
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(pc.sum[2,1:6]*100, 1), "%)", sep = "")
pca.comp$sample <- colnames(vcf.SNPs@gt)[-1]
pca.comp <- merge(pca.comp, samples_data[, -(2:6)], by.x = "sample", by.y="ID")

pca12.plot <- ggplot(pca.comp) +
  geom_point(aes(pca1, pca2, col = Waterbody, shape = Ecotype)) +
  labs(x = pca.labs[1], y = pca.labs[2])
pca23.plot <- ggplot(pca.comp) +
  geom_point(aes(pca2, pca3, col = Waterbody, shape = Ecotype)) +
  labs(x = pca.labs[2], y= pca.labs[3])
pca45.plot <- ggplot(pca.comp) +
  geom_point(aes(pca4, pca5, col = Waterbody, shape = Ecotype)) +
  labs(x = pca.labs[4], y= pca.labs[5])
pca56.plot <- ggplot(pca.comp) +
  geom_point(aes(pca5, pca6, col = Waterbody, shape = Ecotype)) +
  labs(x = pca.labs[5], y= pca.labs[6])

pca.all.plot <- (pca12.plot + pca23.plot)/(pca45.plot + pca56.plot)

ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_paired_PCA.pdf"), pca.all.plot, width = 15, height = 15)
ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_paired_PCA.png"), pca.all.plot, width = 15, height = 15)

# # # # # # # # # # # # # # # #
####### Kinship analysis ######
# # # # # # # # # # # # # # # #

#install.packages("popkin")
library(popkin)

# Load geno file and change missing values to "NA
geno.kin <- geno
geno.kin[geno.kin=="9"] <- NA

vcf.SNP

# Choose subpopulations
#Check dimensions are correct
match(col.names(vcf.SNPs@gt)[-1], samples_data$ID)
paired_samples <- samples_data[match(colnames(vcf.SNPs@gt)[-1], samples_data$ID),]

dim(paired_samples)[1]==dim(geno)[1]
subpops <- paired_samples$Population
# Reorder samples
kin.plot.order <- order(paste(paired_samples$Ecotype ,paired_samples$Population))
paired_samples <- paired_samples[kin.plot.order,]
geno.kin <- geno.kin[kin.plot.order,]
subpops <- subpops[kin.plot.order]
subpops.site.sub <- paired_samples$Population[kin.plot.order]

kinship <- popkin(t(geno.kin), subpops = subpops)
png("test.png", width = 2000, height = 2000)
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
kinship.df$Population_1 <- paired_samples$Population[match(kinship.df$ID.1, paired_samples$ID)]
kinship.df$Population_2 <- paired_samples$Population[match(kinship.df$ID.2, paired_samples$ID)]

kinship.df$Kinship[kinship.df$ID.1==kinship.df$ID.2] <- NA
## Plot
kin.plot <- ggplot(kinship.df) +
  geom_tile(aes(paste(ID.1), paste(ID.2), fill = Kinship)) +
  geom_point(aes(ID.1, 0, col = Population_1), size = 5) +
  geom_point(aes(0, ID.2, col = Population_2), size = 5) +
  scale_fill_gradient2(low = "white", high = "red")
  # theme(axis.text = element_blank())


ggsave("test.png", kin.plot, width = 10, height = 10)
#Calculates structure for samples from K=1 to k=15
max.K <- 6
# MAY NEED TO PAUSE ONEDRIVE
# File names are becoming too Long
obj.at <- snmf(paste0(plot.dir,SNP.library.name,"_paired.geno"), K = 1:max.K, ploidy = 2, entropy = T,
              CPU = 8, project = "new", repetitions = 10, alpha = 100)
stickleback.snmf <- load.snmfProject(file = paste0(plot.dir,SNP.library.name,"_paired.snmfProject"))
stickleback.snmf.sum <- summary(stickleback.snmf)

plot(stickleback.snmf, col = "blue4", cex = 1.4, pch = 19)

ce <- cbind(1:max.K, t(stickleback.snmf.sum$crossEntropy))
colnames(ce) <- c("K", "min","mean","max")
ce <- data.frame(ce)

summary(stickleback.snmf)

which.min(ce$mean)
#Choose K
K <- which.min(ce$mean)
# K <- 4
ce.plot <- ggplot(ce) +
  geom_point(aes(K, mean), size = 2, shape = 19) +
  geom_errorbar(aes(x = K, ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +
  ggtitle(paste("Minimum mean cross-entropy value K =", which.min(ce$mean))) +
  ylab("Cross-entropy") +
  theme_bw()
ce.plot
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_paired_LEA_K1-",max.K,"_cross_entropy.pdf"), plot = ce.plot)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_paired_LEA_K1-",max.K,"_cross_entropy.jpg"), plot = ce.plot)

best <- which.min(cross.entropy(stickleback.snmf, K = K))
qmatrix = Q(stickleback.snmf, K = K, run = best)
dim(qmatrix)
# Tidy data for plotting
qtable <-  cbind(colnames(vcf.SNPs@gt)[-1], rep(1:K, each = dim(qmatrix)[1]), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","Qid", "Q")

#qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$Country_Ocean, sites$Lat))])
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
  facet_grid(~Population+Ecotype, drop = T, scales = "free", space = "free") +
  ylab(label = paste("K =", K))
v

ggsave(filename = "test.pdf", v, width = 18, height = 6)

ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,".pdf"), v, width = 18, height = 6)
ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,".png"), v, width = 18, height = 6)

pop <- unique(samples_data$Population)

pop
#Number of unique sites
Npop = length(pop)
Npop
# Creating qmatrix
qpop = matrix(NA, ncol = K, nrow = Npop)
qpop
coord.pop = matrix(NA, ncol = 2, nrow = Npop)
for (i in 1:length(unique(pop))){
  tmp.pop <- unique(pop)[i]
  print(samples_data$Population[which(samples_data$Population==tmp.pop)])
  print(paste("There are ", length(which(samples_data$Population==tmp.pop)), "samples"))
  if(length(which(samples_data$Population==tmp.pop)) == 1) {
    qpop[i,] <- qmatrix[samples_data$Population == tmp.pop,]
    coord.pop[i,] <- apply(samples_data[samples_data$Population == tmp.pop,][,c("Lat","Long")], 2, mean)
  } else {
    qpop[i,] = apply(qmatrix[samples_data$Population == tmp.pop,], 2, mean)
    coord.pop[i,] = apply(samples_data[samples_data$Population == tmp.pop,][,c("Lat","Long")], 2, mean)
  }
}

print("Check point 1")
q.coord.pop <- data.frame(pop, coord.pop, qpop)

q.coord.pop

print("Check point 2")

colnames(q.coord.pop) <- c("site", "Lat", "Lon", LETTERS[1:K])

print("Check point 3")
q.coord.pop$Lat <- as.numeric(q.coord.pop$Lat)
print("Check point 4")
q.coord.pop$Lon <- as.numeric(q.coord.pop$Lon)
print("Check point 5")


######### 
######### 
######### pca.data <- cbind(sites, pca.comp)
######### #Random order for ploting
######### pca.data <- pca.data[sample(1:length(pca.data$samples)),]
######### #plot(pca.comp)
######### 
######### 
######### ########
######### #PLOTS##
######### ########
######### 
######### world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
######### mex.e <- data.frame( Long = c(-105,-88), Lat = c(23,15))
######### mex.e.zoom <- data.frame( Long = c(-96,-94), Lat = c(18,16))
######### CR.e <- data.frame(Long = c(-85,-82), Lat = c(8,11))
######### 
######### 
######### source("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/theme_black.R")
######### source("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/World/hydrobasins_extract_code.R")
######### 
######### # Colour blind pallette
######### # cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")
######### 
######### p <- ggplot() +
#########   geom_polygon(data = worldmap, aes(long, lat, group = group), fill = "#FFE6D4", col = "black") +
#########   #geom_point(data = sites, aes(x = Long, y = Lat)) +
#########   #geom_point(data = q.coord.pop, aes(x = Long, y = Lat)) 
#########   geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, group = site, r = 1), cols = LETTERS[1:K]) +
#########   #theme_bw()  +
#########   #xlim(c(-130,-60)) +
#########   #ylim(c(-0,60)) +
#########   #coord_sf(expand = FALSE) +
#########   scale_fill_manual(values = cbPalette) +
#########   theme(strip.text = element_text(face = "italic"),
#########         legend.text = element_text(face = "italic")) +
#########   labs(x = "Longitude", y = "Latitude", col = "Species") +
#########   theme(panel.background = element_rect(fill = NA),
#########         panel.ontop = TRUE) +
#########   theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
#########         panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
#########   theme(legend.position = "none") +
#########   coord_map("bonne", Lat0 = 50, xlim = world.e$Long, ylim = world.e$Lat) +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm"))  +
#########   ggtitle(paste("(a) K = ", K))
######### 
######### q <- ggplot() +
#########   geom_raster(data = hill.df.Mex, aes(lon, lat, fill = hill), alpha = 1) +
#########   #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
#########   scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
#########   new_scale_fill() +
#########   geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
#########   geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
#########   geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.2 , group = site), cols = LETTERS[1:K]) +
#########   theme(legend.position="none") +
#########   #theme_bw()  +
#########   #xlim(c(-130,-60)) +
#########   #ylim(c(-0,60)) +
#########   #coord_sf(expand = FALSE) +
#########   scale_fill_manual(values = cbPalette) +
#########   theme(strip.text = element_text(face = "italic"),
#########         legend.text = element_text(face = "italic")) +
#########   labs(x = "Longitude", y = "Latitude", col = "Species") +
#########   theme(panel.background = element_rect(fill = NA),
#########         panel.ontop = TRUE) +
#########   theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
#########         panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
#########   theme(legend.position = "none") +
#########   coord_sf(xlim = mex.e$Long, ylim = mex.e$Lat) +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
#########   ggtitle("(b)")
######### 
######### q.zoom <- ggplot() +
#########   geom_raster(data = hill.df.Mex.z, aes(lon, lat, fill = hill), alpha = 1) +
#########   #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
#########   scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
#########   new_scale_fill() +
#########   geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
#########   geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
#########   geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.08 , group = site), cols = LETTERS[1:K]) +
#########   theme(legend.position="none") +
#########   #theme_bw()  +
#########   #xlim(c(-130,-60)) +
#########   #ylim(c(-0,60)) +
#########   #coord_sf(expand = FALSE) +
#########   scale_fill_manual(values = cbPalette) +
#########   theme(strip.text = element_text(face = "italic"),
#########         legend.text = element_text(face = "italic")) +
#########   labs(x = "Longitude", y = "Latitude", col = "Species") +
#########   theme(panel.background = element_rect(fill = NA),
#########         panel.ontop = TRUE) +
#########   theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
#########         panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
#########   theme(legend.position = "none") +
#########   coord_sf(xlim = mex.e.zoom$Long, ylim = mex.e.zoom$Lat) +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
#########   ggtitle("(b)")
######### 
######### r <- ggplot() +
#########   geom_raster(data = hill.df.CR, aes(lon, lat, fill = hill), alpha = 1) +
#########   scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
#########   theme(legend.position = "none") +
#########   new_scale_fill() +
#########   #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
#########   #scale_fill_gradientn(colours=c("#d95f0e","#FFE6D4")) +
#########   new_scale_color() +
#########   geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
#########   geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
#########   #geom_point(data = sites, aes(x = Long, y = Lat)) +
#########   #geom_point(data = q.coord.pop, aes(x = Long, y = Lat)) 
#########   geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.2 , group = site), cols = LETTERS[1:K]) +
#########   theme(legend.position="none") +
#########   #theme_bw()  +
#########   #xlim(c(-130,-60)) +
#########   #ylim(c(-0,60)) +
#########   #coord_sf(expand = FALSE) +
#########   scale_fill_manual(values = cbPalette) +
#########   theme(strip.text = element_text(face = "italic"),
#########         legend.text = element_text(face = "italic")) +
#########   labs(x = "Longitude", y = "Latitude", col = "Species") +
#########   theme(panel.background = element_rect(fill = NA),
#########         panel.ontop = TRUE) +
#########   theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
#########         panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
#########   theme(legend.position = "none") +
#########   coord_sf(xlim = CR.e$Long, ylim = CR.e$Lat) +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
#########   ggtitle("(c)")
######### 
######### best <- which.min(cross.entropy(stickleback.snmf, K = K))
######### qmatrix = Q(stickleback.snmf, K = K, run = best)
######### # qtable <- cbind(rep(sites$samples,K), rep(sites$Lat,K), rep(1:K, each = length(sites$samples)), c(qmatrix[,1:K]))
######### qtable <-  cbind(colnames(vcf.SNPs@gt)[-1], rep(1:K, each = dim(qmatrix)[1]), c(qmatrix[,1:K]))
######### qtable <-  data.frame(qtable)
######### colnames(qtable) <- c("sample","Qid", "Q")
######### 
######### #qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$Country_Ocean, sites$Lat))])
######### qtable$Q <- as.numeric(qtable$Q)
######### #qtable$sample_Lat <- paste(qtable$Lat, qtable$sample, sep = "_")
######### 
######### #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")
######### cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")
######### 
######### v <- ggplot(qtable)+
#########   geom_bar(stat="identity", aes(sample, Q, fill = as.factor(Qid)), position = "stack", width = 1, col = "black") +
#########   scale_fill_manual(values = cbPalette) +
#########   theme_classic() +
#########   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#########   theme(legend.position = "none") +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
#########   ylab(label = paste("K =", K))+
#########   ggtitle("(d)")
######### v
######### ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"_snp",snp_sub_text,".pdf"), v, width = 18, height = 6)
######### ggsave(filename = paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"_snp",snp_sub_text,".png"), v, width = 18, height = 6)
######### 
######### #Creates multiple K barcharts
######### 
######### #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")
######### 
######### s <- list()
######### 
######### max.K <- min(c(max.K, 6))
######### 
######### for(i in 2:max.K){
#########   best <- which.min(cross.entropy(stickleback.snmf, K = i))
#########   qmatrix = Q(stickleback.snmf, K = i, run = best)
#########   qtable <- cbind(rep(sites$samples,i), rep(sites$Lat,i), rep(1:i, each = length(sites$samples)), c(qmatrix[,1:i]))
#########   qtable <-  data.frame(qtable)
#########   colnames(qtable) <- c("sample","Lat","Qid", "Q")
#########   
#########   qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$Country_Ocean, sites$Lat))])
#########   qtable$Q <- as.numeric(qtable$Q)
#########   qtable$sample_Lat <- paste(qtable$Lat, qtable$sample, "_")
#########   
#########   s[[i]] <- ggplot(qtable)+
#########     geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
#########     scale_fill_manual(values = cbPalette) +
#########     theme_classic() +
#########     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#########     theme(legend.position = "none") +
#########     theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
#########     theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
#########     ylab(label = paste("K =",i))
#########   
######### }
######### 
######### 
######### s[[max.K]] <- ggplot(qtable)+
#########   geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
#########   scale_fill_manual(values = cbPalette) +
#########   theme_classic() +
#########   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#########   theme(legend.position = "none") +
#########   theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
#########   ylab(label = paste("K =",max.K))
######### 
######### 
######### plot2 <- s[[2]]
######### for(i in 3:max.K){
#########   plot2 <- plot2 / s[[i]]
######### }
######### s <- plot2
######### 
######### ## PCA plots
######### 
######### # cbPalette <- c("#999999", "#009E73", "#0072B2", "#E69F00", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7", "black")
######### plot(1:length(cbPalette), col = cbPalette, pch = 19, cex =8)
######### x <- ggplot(pca.data) +
#########   #geom_point(aes(pca1, pca2, colour = Country_Ocean), size = 4) +
#########   #scale_colour_manual(values = cbPalette, name = "Region") +
#########   geom_text(aes(pca1, pca2, colour = species.country.drainage, label = samples), alpha = .8) +
#########   #geom_label_repel(aes(pca1, pca2, label = samples, colour = species), size = 0.5
#########   #                 ,box.padding = 0.1,label.padding = 0.1, min.segment.length = 10) +
#########   xlab(pca.labs[1]) +
#########   ylab(pca.labs[2]) + 
#########   theme_bw() 
######### 
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_",snp_sub_text,"_pca1_2.pdf"), plot=x, height=10, width=10)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_",snp_sub_text,"_pca1_2.jpg"), plot=x, height=10, width=10)
######### 
######### #theme(legend.position = "none")
######### y <- ggplot(pca.data) +
#########   geom_point(aes(pca1, pca2, colour = species.country.drainage), size = 6, alpha = 0.5) +
#########   scale_colour_manual(values = cbPalette, name = "Region") +
#########   xlab(pca.labs[1]) +
#########   ylab(pca.labs[2]) +
#########   theme(legend.position = c(0.21, 0.8)) 
######### z <- ggplot(pca.data) +
#########   geom_point(aes(pca3, pca4, colour = species.country.drainage), size = 6, alpha = 0.5) + 
#########   scale_colour_manual(values = cbPalette, name = "Region") +
#########   xlab(pca.labs[3]) +
#########   ylab(pca.labs[4]) +
#########   theme(legend.position = "none") 
######### w <- ggplot(pca.data) +
#########   geom_point(aes(pca5, pca6, colour = species.country.drainage), size = 6, alpha = 0.5) + 
#########   scale_colour_manual(values = cbPalette, name = "Region") +
#########   #xlab(pca.labs[5]) +
#########   #ylab(pca.labs[6]) +
#########   theme(legend.position = "none")
######### ww <- ggplot(pca.data) +
#########   geom_label(aes(pca5, pca6, label = samples, colour = species.country.drainage), size = 6, alpha = 0.5) +
#########   scale_colour_manual(values = cbPalette, name = "Region") +
#########   xlab(pca.labs[5]) +
#########   ylab(pca.labs[6]) +
#########   theme(legend.position = "none")
######### a <- ggplot() +
#########   geom_polygon(data = worldmap, aes(long, lat, group = group), fill = "#FFE6D4", col = "black") +
#########   geom_point(data = pca.data, aes(Long, Lat, colour = species.country.drainage), size = 6) +
#########   scale_colour_manual(values = cbPalette, name = "Region") +
#########   scale_fill_manual(values = cbPalette) +
#########   theme(strip.text = element_text(face = "italic"),
#########         legend.text = element_text(face = "italic")) +
#########   labs(x = "Longitude", y = "Latitude", col = "Species") +
#########   theme(panel.background = element_rect(fill = NA),
#########         panel.ontop = TRUE) +
#########   theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
#########         panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
#########   theme(legend.position = "none") +
#########   coord_map("bonne", Lat0 = 50, xlim = world.e$Long, ylim = world.e$Lat)
######### 
######### 
######### plot1 <- (p+q+r)/v + plot_layout(heights = c(4, 1))
######### plot2 <- y+z + w + a + plot_layout(nrow = 2)
######### 
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.pdf"), plot=plot1, height=7, width=15)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.jpg"), plot=plot1, height=7, width=15)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_Mex_z.pdf"), plot=q.zoom, height=10, width=10)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_Mex_z.jpg"), plot=q.zoom, height=10, width=10)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_barplot_1-",max.K,"_snp",snp_sub_text,"_complete_DP10_basins_patchwork.pdf"), plot=s, height=20, width=15)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_barplot_1-",max.K,"_snp",snp_sub_text,"_complete_DP10_basins_patchwork.jpg"), plot=s, height=20, width=15)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_complete","snp",snp_sub_text,"_DP10_basins_patchwork.pdf"), height=9, width=10, plot=plot2)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_complete","snp",snp_sub_text,"_DP10_basins_patchwork.jpg"), height=9, width=10, plot=plot2)
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_1-2-complete","snp",snp_sub_text,"_DP10_basins_patchwork.pdf"), height=8, width=5, plot=(y/a))
######### ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_1-2-complete","snp",snp_sub_text,"_DP10_basins_patchwork.jpg"), height=8, width=5, plot=(y/a))