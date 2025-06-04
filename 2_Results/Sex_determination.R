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

sex_chr <- c("NC_053230.1","NC_053233.1")

chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T))
select_chr <- c("NC_053212.1","NC_053230.1","NC_053233.1")

select_chr <- chr[chr$RefSeq.seq.accession%in%select_chr,c("Sequence.name","RefSeq.seq.accession")]
print(chr[,c("Sequence.name","RefSeq.seq.accession")], n = 22)

sex_chr <- c("NC_053230.1","NC_053233.1")

vcf.SNPs <- list()
## Loop through and convert raw bcfs into vcf for each 
for(i in 1:length(select_chr$Sequence.name)){
    print(select_chr$Sequence.name[i])
    if(!paste0("stickleback_",select_chr$RefSeq.seq.accession[i], "_sorted.vcf.gz")%in%list.files(paste0(dir.path, "/vcfs/"))){
      print("vcf does not exist: converting from bcf to vcf")
      system(paste0("bcftools view -O z ",dir.path, "/vcfs/stickleback_", select_chr$RefSeq.seq.accession[i], "_sorted.bcf > ",dir.path, "/vcfs/stickleback_",select_chr$RefSeq.seq.accession[i], "_sorted.vcf.gz"))
    }
    vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]] <- read.vcfR(paste0(dir.path, "/vcfs/stickleback_",select_chr$RefSeq.seq.accession[i], "_sorted.vcf.gz"),
                      verbose = T)
    vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]] <- (vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]])[samples = sort(colnames(vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]])[-1])] 

}


# Get an read sample information
samples_data <- data.frame(ID = colnames(vcf.SNPs[[1]]@gt)[-1])
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)
samples_data <- samples_data[samples_data$ID%in%colnames(vcf.SNPs[[1]]@gt),]
samples_data <- samples_data[match(samples_data$ID, (colnames(vcf.SNPs[[1]]@gt)[-1])),]

myDepth <- lapply(vcf.SNPs, function(x) {
extract.gt(x, element ="DP")
})

## Cbind all scaffolds depths into single data frame
myDepth.df <- lapply(myDepth, function(x) {
    scaf.tmp <- paste0("NC_", stringr::str_split_i(rownames(x), "_", c(2)))
    pos.tmp <- stringr::str_split_i(rownames(x), "_", c(3))
return(cbind.data.frame(scaf = scaf.tmp, pos = pos.tmp, x))
})

# Convert to data frame
myDepth.df <- do.call("rbind", myDepth.df)

## Reformat columns
myDepth.df.long <- pivot_longer(myDepth.df, cols = colnames(vcf.SNPs[[1]]@gt)[-1], values_to = "depth", names_to = "ID")
myDepth.df.long$pos <- as.numeric(myDepth.df.long$pos)
myDepth.df.long$ID <- as.factor(myDepth.df.long$ID)
myDepth.df.long$depth <- as.numeric(myDepth.df.long$depth)

## Calculate mean coverage for each chromosome was each sample
chr_cov <- myDepth.df.long %>%
group_by(ID) %>%
filter(depth<mean(depth)*10) %>%
summarise(avg.chr = mean(depth[!scaf%in%sex_chr],na.rm = T), avg.x = mean(depth[scaf==sex_chr[1]],na.rm = T), avg.y = mean(depth[scaf==sex_chr[2]],na.rm =T))
## Calculate sample ratio of X to Y coverage for each sample
chr_cov$x.ratio <- chr_cov$avg.x/chr_cov$avg.chr
chr_cov$y.ratio <- chr_cov$avg.y/chr_cov$avg.chr

chr_cov <- merge(chr_cov, samples_data, by = "ID")

chr_cov$Sex[is.na(chr_cov$Sex)] <- "?"

## Determin ratio value of X and Y coverage
Xcov_sex_determine <- 0.8
## Create sex column
chr_cov$Gsex <- c("M","F")[as.numeric(chr_cov$x.ratio>=Xcov_sex_determine)+1]
chr_cov$GsexN <- as.numeric(chr_cov$x.ratio>=Xcov_sex_determine)+1

# Save as data frame
write.table(chr_cov, file = paste0(dir.path, "/vcfs/Genomic_sex_determination.csv"), 
            row.names = F, quote = F, col.names = TRUE, sep = ",")

# Save as input for bcftools call function
write.table(chr_cov[,c("ID", "Gsex")], file = paste0("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/Genomic_sex_determination.txt"), 
            row.names = F, quote = F, col.names = FALSE, sep = "\t")

# Save as input for bcftools call function (PED) (Fam1 is ignored by bcftools)
write.table(cbind.data.frame("Fam1", chr_cov[,c("ID")], "0", "0", chr_cov[,c("Gsex")]), file = paste0("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/Genomic_sex_determination.ped"), 
            row.names = F, quote = F, col.names = FALSE, sep = "\t")

p <- ggplot(chr_cov, aes(Sex, x.ratio)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(aes(col = Ecotype), height = 0, width = 0.2) +
geom_hline(yintercept = Xcov_sex_determine, col = "red") +
ylab("Mean X coverage / Mean ChrI coverage")

q <- ggplot(chr_cov, aes(Sex, y.ratio)) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(aes(col = Ecotype), height = 0, width = 0.2) +
ylab("Mean Y coverage / Mean ChrI coverage")

z <- ggplot(chr_cov, aes(x.ratio)) +
geom_histogram(fill = "lightblue") +
geom_vline(xintercept = Xcov_sex_determine, col = "red") +
xlab("mean X coverage / Mean ChrI coverage")

ggsave(paste0(plot.dir,"Genomic_sex_determination.pdf"), (p+q)/z + plot_layout(guides = "collect"), height= 7, width = 12)
ggsave(paste0(plot.dir,"Genomic_sex_determination.png"), (p+q)/z + plot_layout(guides = "collect"), height= 7, width = 12)

## Rolling coverage across the genome between sexs
library(zoo)

## Add sex status to Depth dataframe
myDepth.df.long$Gsex <- chr_cov$Gsex[match(myDepth.df.long$ID, chr_cov$ID)]
myDepth.df.long$avg.chr <- chr_cov$avg.chr[match(myDepth.df.long$ID, chr_cov$ID)]
myDepth.df.long$avg.chr.ratio <- myDepth.df.long$depth/myDepth.df.long$avg.chr

## Set bin length
bin.length <- 250000
# List to add info too
bin.Mean <- list()
for(chr.i in select_chr$RefSeq.seq.accession){
  start <- 1
  while(start < chr$Seq.length[chr$RefSeq.seq.accession==chr.i]){
    end <- start + bin.length
    win.tmp <- paste0(start,"-", end, "-",chr.i)
    print(paste0(win.tmp))
    ## Calculated mean
    # subset to window
    win.vcf.tmp <- myDepth.df.long[myDepth.df.long$pos>=start & myDepth.df.long$pos<end & myDepth.df.long$scaf==chr.i,]
    # Calculate mean coverage of M and F
    bin.Mean[[win.tmp]] <- c(F.mn = mean(win.vcf.tmp$avg.chr.ratio[win.vcf.tmp$Gsex=="F"]), M.mn = mean(win.vcf.tmp$avg.chr.ratio[win.vcf.tmp$Gsex=="M"]))
    # Shift along window
    start <- start + bin.length
    
  }
}

# Convert to data frame
bin.Mean.df <- as.data.frame(do.call("rbind",bin.Mean))
## Add in window size
bin.Mean.df <- cbind(bin.Mean.df, stringr::str_split(rownames(bin.Mean.df), pattern = "-", simplify = TRUE))
colnames(bin.Mean.df) <- c(colnames(bin.Mean.df)[1:2], "start", "end", "scaf")
bin.Mean.df$F.mn <- as.numeric(bin.Mean.df$F.mn)
bin.Mean.df$M.mn <- as.numeric(bin.Mean.df$M.mn)
bin.Mean.df$start <- as.numeric(bin.Mean.df$start)

## Add chr name
bin.Mean.df$chr <- chr$Chromosome.name[match(bin.Mean.df$scaf, chr$RefSeq.seq.accession)]

p <- ggplot(bin.Mean.df) +
    geom_hline(yintercept = c(0.5,1), linetype = "dashed") +
    facet_grid(.~chr) +
    ylim(c(0,2)) +
    theme_bw()

ggsave(paste0(plot.dir,"Genomic_sex_depth_window_ylim.pdf"), p, width = 10, height = 7)
ggsave(paste0(plot.dir,"Genomic_sex_depth_window_ylim.png"), p, width = 10, height = 7)

p <- ggplot(bin.Mean.df) +
    geom_hline(yintercept = c(0.5,1), linetype = "dashed") +
    geom_point(aes(start, F.mn, col = "Genomic Female")) +
    geom_point(aes(start, M.mn,col = "Genomic Male")) +
    scale_color_manual(values = c("orange", "skyblue")) +
    facet_grid(.~chr) +
    theme_bw() +
    theme(legend.position = "top")

ggsave(paste0(plot.dir,"Genomic_sex_depth_window.pdf"), p, width = 10, height = 7)
ggsave(paste0(plot.dir,"Genomic_sex_depth_window.png"), p, width = 10, height = 7)

q <- ggplot(bin.Mean.df) +
    geom_hline(yintercept = c(0, 0.5, -0.5), linetype = "dashed") +
    geom_point(aes(start, F.mn-M.mn)) +
    facet_grid(.~chr) +
    theme_bw()


ggsave(paste0(plot.dir,"Genomic_sex_depth_window.pdf"), p/q, width = 10, height = 10)
ggsave(paste0(plot.dir,"Genomic_sex_depth_window.png"), p/q, width = 10, height = 10)
ggsave(paste0("test.pdf"), p/q, width = 10, height = 10)
