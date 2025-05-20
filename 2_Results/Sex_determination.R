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
    # system(paste0("bcftools view -O z ",dir.path, "/vcfs/stickleback_", select_chr$RefSeq.seq.accession[i], "_sorted.bcf > ",dir.path, "/vcfs/stickleback_",select_chr$RefSeq.seq.accession[i], "_sorted.vcf.gz"))
    vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]] <- read.vcfR(paste0(dir.path, "/vcfs/stickleback_",select_chr$RefSeq.seq.accession[i], "_sorted.vcf.gz"),
                      verbose = T)
    vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]] <- (vcf.SNPs[[select_chr$RefSeq.seq.accession[i]]])[samples = sort(colnames(vcf.SNPs.Y@gt)[-1])] 

}


# Get an read sample information
samples_data <- data.frame(ID = colnames(vcf.SNPs[[1]]@gt)[-1])
samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
samples_data <- merge(samples_data, samples, by.x = "ID", by.y="individual", all.x = T)
samples_data <- samples_data[samples_data$ID%in%colnames(vcf.SNPs[[1]]@gt),]
samples_data <- samples_data[match(samples_data$ID, (colnames(vcf.SNPs[[1]]@gt)[-1])),]
print("Do any samples names not line up with (False is good)")
any(!samples_data$ID==(colnames(vcf.SNPs.X@gt)[-1]))

myDepth <- lapply(vcf.SNPs, function(x) {
extract.gt(x, element ="DP")
})

myDepth.df <- lapply(myDepth, function(x) {
    scaf.tmp <- paste0("NC_", stringr::str_split_i(rownames(x), "_", c(2)))
    pos.tmp <- stringr::str_split_i(rownames(x), "_", c(3))
return(cbind.data.frame(scaf = scaf.tmp, pos = pos.tmp, x))
})

myDepth.df <- do.call("rbind", myDepth.df)
 
myDepth.df.long <- pivot_longer(myDepth.df, cols = colnames(vcf.SNPs[[1]]@gt)[-1], values_to = "depth", names_to = "ID")
myDepth.df.long$pos <- as.numeric(myDepth.df.long$pos)
myDepth.df.long$ID <- as.factor(myDepth.df.long$ID)
myDepth.df.long$depth <- as.numeric(myDepth.df.long$depth)

chr_cov <- myDepth.df.long %>%
group_by(ID) %>%
filter(depth<mean(depth)*10) %>%
summarise(avg.chr = mean(depth[!scaf%in%sex_chr],na.rm = T), avg.x = mean(depth[scaf==sex_chr[1]],na.rm = T), avg.y = mean(depth[scaf==sex_chr[2]],na.rm =T))
chr_cov$x.ratio <- chr_cov$avg.x/chr_cov$avg.chr
chr_cov$y.ratio <- chr_cov$avg.y/chr_cov$avg.chr

chr_cov <- merge(chr_cov, samples_data, by = "ID")

chr_cov$Sex[is.na(chr_cov$Sex)] <- "?"

Xcov_sex_determine <- 0.8

chr_cov$Gsex <- c("M","F")[as.numieric(chr_cov$x.ratio>=Xcov_sex_determine)+1]

write.table(chr_cov, file = paste0(dir.path, "/vcfs/Genomic_sex_determination.csv"), 
            row.names = F, quote = F, col.names = TRUE, sep = ",")

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


