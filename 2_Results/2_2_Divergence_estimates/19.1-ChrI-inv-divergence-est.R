library(vcfR)
library(LEA)
library(data.table)
library(tidyverse)

# Read in vcf
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
# vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz"
vcf.SNPs <- read.vcfR(vcf.file)
outdir <- dirname(vcf.file)

## Get coding region of atp1a1a
# Read vcf
vcf.atp1a1a.file <- args[2]
# vcf.atp1a1a.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.vcf.gz"
vcf.atp1a1a <- read.vcfR(vcf.atp1a1a.file)

# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[, c("FORMAT", sort(colnames(vcf.SNPs@gt)[-1]))] 

## Get list of chromosomes
chr <- unique(vcf.SNPs@fix[,"CHROM"])

# Read in atp1a1a
atp1a1a.bed <- read.table(paste0(outdir, "/atp1a1a_CDS.bed")) %>%
    rename(chr = V1, start = V2, end = V3)

# Set window size
wndsize <- 25000
wndslid <- 5000

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

      # Calculate the number of heterozygousity
      het.count <- apply(geno.mat, MARGIN = 2, function(x) sum(x[(!is.na(x))]=="1"))[colnames(geno.mat)]

      # Check is dist matrix contains any NA
      if(any(is.na(dc))) {return(NULL) ; print("mds-failed")}
      mds <- cmdscale(dc, k = 2)

      # Return data.frame with PCA1, PCA2, MDS1, MDS2, sample names, window name
      return(data.frame(sample = colnames(geno.mat), chr = chr, start = wnd$start, end = wnd$end, base.calls = base.calls, nsamps = ncol(geno.mat),
                        MDS1 = mds[,1], MDS2 = mds[,2], dist.marine = dc.marine.samp, het.prop = het.count/base.calls))

    }
  )


samples <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
# Override Ecotypes
samples$Ecotype[samples$Ecotype=="st"] <- "fw"
samples$Ecotype[samples$Ecotype=="anad"] <- "mig"

# Merge mds with sample data
samples_data <- merge(mds.comp, samples,  by.y = "individual", by.x="sample", all.x = T) %>%
                  mutate(windowname = paste(chr, start, end, sep = "-"),
                         mid = end - ((end-start)/2))

# Transform PCA so that the axis is also segregating populations in the same direction across all windows
## Copy over PCA and MDS data to new scaled columns
samples_data$MDS1_scaled <- samples_data$MDS1
samples_data$MDS2_scaled <- samples_data$MDS2

# Define the columns you want to check and potentially invert
scale_cols <- c("MDS1_scaled", "MDS2_scaled")

# Compute the sign for each window and ecotype == "mig"
signs <- samples_data %>%
  filter(Ecotype == "mig") %>%
  group_by(windowname) %>%
  summarise(across(all_of(scale_cols), ~ sign(median(.x, na.rm = TRUE)), .names = "sign_{.col}"), .groups = "drop")

# Join sign info back to original data
samples_data <- samples_data %>%
  left_join(signs, by = "windowname") %>%
  mutate(across(all_of(scale_cols),
                ~ ifelse(get(paste0("sign_", cur_column())) == -1, -.x, .x))) %>%
  select(-starts_with("sign_"))  # remove helper columns

# Calculate max and min MDS
max_min <- samples_data %>%
  group_by(windowname) %>%
  summarise(
    max.mds = abs(max(MDS1_scaled, na.rm = TRUE)),
    min.mds = abs(min(MDS1_scaled, na.rm = TRUE)),
    .groups = 'drop'
  )

## Calculate position of individual samples within max and min MDS
samples_data <- samples_data %>%
  left_join(max_min, by = "windowname") %>%
  mutate(
    MDS1_ratio = (MDS1_scaled + min.mds) / (max.mds + min.mds),
    Ecotype <- factor(Ecotype, levels = c("mig", "resi", "fw")) 
  )


min(c(atp1a1a.bed$start,atp1a1a.bed$end))
max(c(atp1a1a.bed$start,atp1a1a.bed$end))


mds_chrI_inv <- ggplot(samples_data) +
  geom_tile(aes(as.numeric(mid), sample, fill = MDS1_ratio)) +
  # geom_vline(data = atp1a1a.bed, aes(xintercept = end), alpha = 0.5) +
  geom_vline(xintercept = as.numeric(c(min(atp1a1a.bed$start), max(atp1a1a.bed$end))), col = "black") +
  # geom_segment(data = atp1a1a.bed, aes(x = start, xend = end, y = 1, col = transcript.id), linewidth = 4, alpha = 0.5) +
  scale_fill_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", midpoint=0.5, name =  "MDS Scaled") +
  # scale_fill_gradient2(low = "#009E73", mid = "#E69F00", high = "#56B4E9", midpoint=0.5, name =  "MDS Scaled") +
  #scale_fill_gradient2(low = "firebrick3", mid = "orange" ,high = "darkgreen", midpoint=0.5, name =  "MDS Scaled") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6),name = "Mbps", expand = c(0.01,0)) +
  facet_grid(Ecotype+Population~.,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "bottom", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(0.5,'lines'),
        legend.frame = element_rect(colour="black"),
        legend.ticks = element_line(colour="black"),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
        # axis.text.y = element_blank(),                    # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        axis.title.y = element_blank(), 
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(fill = "grey90", color = "black", size = 0.5))

ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"),
       plot = mds_chrI_inv, width = 30, height = 30)

ggsave(paste0("test.png"),
       plot = mds_chrI_inv, width = 30, height = 30)


p <- ggplot(samples_data) +
    geom_line(aes(mid, dist.marine/base.calls, col = Population, group = sample)) +
    geom_vline(xintercept = c(26836909, 26867066)) +
    annotate("text", x = (26836909+((26867066-26836909)/2)), y = 0.0, label = "atp1a1a", color = "black", size = 4, vjust = 1) +
    geom_segment(data = atp1a1a.bed, aes(x = start, xend = end, y = -0.0001), linewidth = 10) +
    theme_bw() +
    coord_cartesian(clip = "off")  # allow drawing in the margin area

p <- ggplot(samples_data) +
    geom_line(aes(mid, MDS1_ratio, col = het.prop, group = sample)) +
    geom_vline(xintercept = c(26836909, 26867066)) +
    annotate("text", x = (26836909+((26867066-26836909)/2)), y = 0.0, label = "atp1a1a", color = "black", size = 4, vjust = 1) +
    scale_color_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", name =  "Heterozygousity") +
    geom_segment(data = atp1a1a.bed, aes(x = start, xend = end, y = -0.0001), linewidth = 10) +
    theme_bw() +
    coord_cartesian(clip = "off") 

## ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"), p, width = 20)
ggsave(paste0("test.png"), p, width = 20)

# Extract atp1a1a genotype
geno.mat.atp1a1a <- extract.gt(vcf.atp1a1a,element = "GT")

# Remove non Uist samples
geno.mat.atp1a1a <- geno.mat.atp1a1a[,colnames(geno.mat.atp1a1a)%in%(samples$individual[samples$Region=="Uist"])]

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
# Heterozygousity
het.count <- apply(geno.mat.atp1a1a, MARGIN = 2, function(x) sum(x[(!is.na(x))]=="1"))[colnames(geno.mat.atp1a1a)]

# Run mds
mds.atp1a1a <- cmdscale(dc.atp1a1a, k = 2)

# Combine
mds.atp1a1a.df <- data.frame(sample = colnames(geno.mat.atp1a1a), base.calls = base.calls.atp1a1a, nsamps = ncol(geno.mat.atp1a1a),
                        MDS1 = mds.atp1a1a[,1], MDS2 = mds.atp1a1a[,2], dist.marine = dc.atp1a1a.marine.samp, het.prop = het.count/base.calls.atp1a1a) %>%
                         merge(samples,  by.y = "individual", by.x="sample", all.x = T)


p.het.atp1a1a <- ggplot(mds.atp1a1a.df) +
  geom_point(aes(MDS1, het.prop, col = Waterbody, shape = Ecotype), size = 3) +
            scale_shape(na.translate = TRUE, na.value = 23)

p.mds.atp1a1a <- ggplot(mds.atp1a1a.df) +
  geom_point(aes(MDS1, MDS2, col = Waterbody, shape = Ecotype), size = 3) +
            scale_shape(na.translate = TRUE, na.value = 23)


#  Calc nj 
library(ape)
library(ggtree)
library(patchwork)
nj.data <- nj(dc.atp1a1a)

# Create tree plot
plot.tree <- ggtree(nj.data, aes(color = Waterbody), layout = "ape")
#  Combine with sample data
plot.tree <- plot.tree %<+% mds.atp1a1a.df
## Custom tip colours
plot.tree <- plot.tree + geom_tippoint(aes(color = Waterbody, shape = Ecotype), stroke = 1, size=3) +
  geom_tiplab(aes(label = paste(Population, label)),hjust = -0.25, size = 2) +
  scale_shape(na.translate = TRUE, na.value = 21) + 
  scale_x_continuous(expand = c(0.25,0.25)) +
  scale_y_continuous(expand = c(0.25,0.25)) +
  ggtitle("Atp1a1a") +
  theme(legend.position = "bottom") 

ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"), 
      mds_chrI_inv + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20)
ggsave(paste0("test.png"), mds_chrI_inv + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20, height = 15)
ggsave(paste0("test.pdf"), mds_chrI_inv + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20, height = 15)


### Write out phyfile
#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes,
# Extract atp1a1a genotype
gt.apt1a1a <- extract.gt(vcf.atp1a1a, element = "GT", return.alleles = T)
# Remove non Uist samples
gt.apt1a1a <- gt.apt1a1a[,colnames(gt.apt1a1a)%in%(samples$individual[samples$Region=="Uist"])]

gt.apt1a1a[gt.apt1a1a=="A/A"] <- "A"
gt.apt1a1a[gt.apt1a1a=="T/T"] <- "T"
gt.apt1a1a[gt.apt1a1a=="G/G"] <- "G"
gt.apt1a1a[gt.apt1a1a=="C/C"] <- "C"
gt.apt1a1a[gt.apt1a1a=="A/G"] <- "R"
gt.apt1a1a[gt.apt1a1a=="G/A"] <- "R"
gt.apt1a1a[gt.apt1a1a=="C/T"] <- "Y"
gt.apt1a1a[gt.apt1a1a=="T/C"] <- "Y"
gt.apt1a1a[gt.apt1a1a=="A/C"] <- "M"
gt.apt1a1a[gt.apt1a1a=="C/A"] <- "M"
gt.apt1a1a[gt.apt1a1a=="G/T"] <- "K"
gt.apt1a1a[gt.apt1a1a=="T/G"] <- "K"
gt.apt1a1a[gt.apt1a1a=="C/G"] <- "S"
gt.apt1a1a[gt.apt1a1a=="G/C"] <- "S"
gt.apt1a1a[gt.apt1a1a=="A/T"] <- "W"
gt.apt1a1a[gt.apt1a1a=="T/A"] <- "W"
gt.apt1a1a[gt.apt1a1a=="."] <- "-"

# Check none of the SNPs are entirely heterozgous and remove them if they are
no.longer.poly <- apply(gt.apt1a1a, MARGIN = 1, function(x) length(unique(x[x!="-"]))>1)
gt.apt1a1a <- gt.apt1a1a[no.longer.poly,]

### ## Write out file
ape::write.dna(t(gt.apt1a1a), file = paste0(outdir,"/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.phy"), format ="interleaved")