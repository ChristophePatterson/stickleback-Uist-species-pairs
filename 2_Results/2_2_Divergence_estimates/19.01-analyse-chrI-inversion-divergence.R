library(vcfR)
library(LEA)
library(data.table)
library(tidyverse)
library(ggnewscale)

# Read in vcf
args <- commandArgs(trailingOnly=T)
#vcf.file <- args[1]
vcf.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz"
vcf.SNPs <- read.vcfR(vcf.file)
outdir <- dirname(vcf.file)

## Get coding region of atp1a1a
# Read vcf
#vcf.atp1a1a.file <- args[2]
vcf.atp1a1a.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest/ChrI_Inv/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.vcf.gz"
vcf.atp1a1a <- read.vcfR(vcf.atp1a1a.file)

# Make vcf be in alphabetical order
vcf.SNPs <- vcf.SNPs[, c("FORMAT", sort(colnames(vcf.SNPs@gt)[-1]))] 

populations <- c("DUIN", "OBSE", "LUIB", "CLAC", 
                 "DUIM", "OBSM", "LUIM", "CLAM", 
                 "OLAV", "TORM")

## REmove unneeded samples
seq_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
# Override Ecotypes
seq_data$Ecotype[seq_data$Ecotype=="st"] <- "fw"
seq_data$Ecotype[seq_data$Ecotype=="anad"] <- "mig"
seq_data$Population[seq_data$Population=="OLST"] <- "OLAV"
seq_data$Population[seq_data$Population=="TOST"] <- "TORM"


seq_data <- seq_data[seq_data$Population%in%populations&!is.na(seq_data$Population),]
seq_data <- seq_data[seq_data$individual %in% colnames(vcf.SNPs@gt)[-1],]
unique(seq_data$Waterbody)
vcf.SNPs <- vcf.SNPs[samples = seq_data$individual]

# Remove invaient and non-biallelic snps
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs, na.omit = T),]


## Get list of chromosomes
chr <- unique(vcf.SNPs@fix[,"CHROM"])

# Read in atp1a1a
atp1a1a.bed <- read.table(paste0(outdir, "/atp1a1a_CDS.bed")) %>%
    rename(chr = V1, start = V2, end = V3)

# Read in all gene bed file
gene.bed <- read.table(paste0(outdir, "/ChrI_inv_CDS.bed")) %>%
  rename(chr = V1, start = V2, end = V3, type = V4, Duke.gene = V5, transcript.id = V6, gene = V7) %>%
  mutate(seg.len = end - start) %>%
  mutate(gene = ifelse(Duke.gene == "g1118", "atp1a1a", gene),
         mid = end-((end-start)/2))

# Set window size
wndsize <- 25000
wndslid <- 5000

# Create windows for each chromosome
## Blank list
chr.df <- as.data.table(vcf.SNPs@fix)[, .(POS = as.numeric(POS)), by = CHROM]
sldwindows.df <- chr.df[, {
  start <- seq(min(POS), max(POS), by = wndslid)
  end <- start + wndsize + 1
  .(start = start, end = end)
}, by = CHROM]

range.plot <- rbind(c(range(c(gene.bed$start, gene.bed$end)),"bed"),
      c(range(chr.df$POS), "extract.vcf"),
      c(range(sldwindows.df$end), "windows"),
      c(range(as.numeric(vcf.SNPs@fix[,"POS"])), "vcf")) %>%
  as.data.frame() %>%
  rename(min.val = V1, max.val = V2, type = V3) %>%
  mutate(min.val = as.numeric(min.val), max.val = as.numeric(max.val)) %>%
  ggplot() +
  geom_segment(aes(x = min.val, xend = max.val, y = 0), col = "red", linewidth = 5) +
  geom_segment(data = sldwindows.df, aes(x = start, xend = end, y = rnorm(n = nrow(sldwindows.df))))
ggsave(file = "test.png", range.plot)

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

# Merge mds with sample data
samples_data <- merge(mds.comp, seq_data,  by.x = "sample", by.y="individual", all.x = T) %>%
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
    Ecotype = factor(Ecotype, levels = c("fw", "mig", "resi", "gene")) 
  )

min(c(atp1a1a.bed$start,atp1a1a.bed$end))
max(c(atp1a1a.bed$start,atp1a1a.bed$end))

# Create intervals to plot 

pack_intervals <- function(dt, start_col, end_col, buffer, out_col = "track") {
  dt <- as.data.table(copy(dt))
  setorderv(dt, c(start_col, end_col))
  
  last_end <- numeric(0)
  track_id <- integer(nrow(dt))
  
  s <- dt[[start_col]]
  e <- dt[[end_col]]
  
  s <- s - buffer
  e <- e + buffer
  
  for (i in seq_len(nrow(dt))) {
    placed <- FALSE
    if (length(last_end) > 0) {
      for (t in seq_along(last_end)) {
        if (s[i] > last_end[t]) {
          track_id[i] <- t
          last_end[t] <- e[i]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      last_end <- c(last_end, e[i])
      track_id[i] <- length(last_end)
    }
  }
  
  dt[, (out_col) := track_id]
  dt
}

gene.bed.track <- pack_intervals(gene.bed[gene.bed$type=="gene",],
                         "start", "end", 100000 ,"track")
gene.bed.track$gene.col <- rep(c("A","B"), length.out = length(gene.bed.track$track))
gene.bed.track <- gene.bed.track %>%
  mutate(Ecotype = factor("gene"),
         Population = "")

gene.bed <- gene.bed %>% 
        left_join(gene.bed.track[,c("Duke.gene", "gene.col", "track")], by = "Duke.gene") %>%
        mutate(Ecotype = factor("gene"),
         Population = "")

# Facet colours
ecotype_cols <- c("mig" = "#1E88E5", "resi" = "#009E73", "fw" = "#FFC107")

mds_chrI_inv_genes <- ggplot(samples_data) +
  geom_tile(aes(as.numeric(mid), sample, fill = MDS1_ratio)) +
  # geom_vline(xintercept = as.numeric(c(min(gene.bed$start), max(gene.bed$end))), col = "black") +
  geom_segment(data = gene.bed.track,aes(x = as.numeric(mid), y = 0, yend = -track*1.5+0.5), show.legend = F) +
  geom_segment(data = gene.bed.track,aes(x = as.numeric(start), xend = end, y = 0, col = gene.col), linewidth = 1, show.legend = F) +
  geom_segment(data = gene.bed[gene.bed$type=="exon",],aes(x = as.numeric(start), xend = end, y = 0, col = gene.col), linewidth = 5, show.legend = F) +
  geom_text(data = gene.bed.track,aes(x = as.numeric(mid), y = -track*1.5, label = gene), size = 4, show.legend = F) +
  # geom_point(data = gene.bed.track, aes(x = start, y = 0, col = gene.col), show.legend = F) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_fill_gradient2(low = "#FFC107", mid = "#D81B60", high = "#1E88E5", midpoint=0.5, name =  "MDS Scaled", 
      guide = guide_colorbar(title.vjust = 0.6, title.hjust = 100,
        barwidth = unit(8, "cm"),  # Adjust the width of the bar
        barheight = unit(1, "cm")   # Adjust the height of the bar
    )) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6),name = "Mbps", expand = c(0.01,0)) +
  facet_grid(Ecotype+Population~.,scale = "free", space = "free", switch = "y") +
  new_scale_fill() +
  geom_tile(aes(x = max(as.numeric(mid))+wndsize, y = sample, fill = Ecotype), width = wndsize, inherit.aes = FALSE, show.legend = F) +
  scale_fill_manual(values = ecotype_cols, name = "Ecotype Group") +
  theme_classic() +
  theme(legend.position = "top", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(0.5,'lines'),
        legend.frame = element_rect(colour="black"),
        legend.ticks = element_line(colour="black"),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
        axis.text.y = element_blank(),                    # hide left/right axis ticks
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        axis.title.y = element_blank(), 
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(fill = "grey90", color = "black", size = 0.5),
        strip.text.y = element_text(size = 12),text = element_text(size = 20))

ggsave(paste0(outdir,"/sliding_window_mds_chrI_inv_wnd",wndsize ,"_sld",wndslid,".png"),
       plot = mds_chrI_inv_genes, width = 15, height = 15)

ggsave(paste0("test.png"),
       plot = mds_chrI_inv_genes, width = 15, height = 15)


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
# ggsave(paste0("test.png"), p, width = 20)

# Extract atp1a1a genotype
geno.mat.atp1a1a <- extract.gt(vcf.atp1a1a,element = "GT")

# Remove non Uist samples
geno.mat.atp1a1a <- geno.mat.atp1a1a[,colnames(geno.mat.atp1a1a)%in%(seq_data$individual[seq_data$Region=="Uist"])]

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
                         merge(seq_data,  by.y = "individual", by.x="sample", all.x = T) %>%
                         mutate(p.miss = base.calls/nrow(geno.mat.atp1a1a))

# Samples with heterozygousity over 0.2
high.het.ind <- mds.atp1a1a.df$sample[mds.atp1a1a.df$het.prop>=0.2]
high.cov.ind <- mds.atp1a1a.df$sample[(mds.atp1a1a.df$p.miss<=0.8)]
print(high.cov.ind)

p.het.atp1a1a <- ggplot(mds.atp1a1a.df) +
  geom_point(aes(MDS1, het.prop, col = Ecotype, shape = Ecotype), size = 3) +
            scale_shape(na.translate = TRUE, na.value = 23) +
  theme_bw()

p.mds.atp1a1a <- ggplot(mds.atp1a1a.df) +
  geom_point(aes(MDS1, MDS2, col = Ecotype, shape = Ecotype), size = 3) +
            scale_shape(na.translate = TRUE, na.value = 23) +
  theme_bw()

atp1a1a.df.geno <- cbind.data.frame(vcf.atp1a1a@fix[,c("CHROM", "POS")], geno.mat.atp1a1a) %>%
  pivot_longer(cols = -c("CHROM", "POS"), values_to = "genotype", names_to = "sample") %>%
  left_join(seq_data[,c("individual", "Waterbody","Population", "Ecotype")], by = c("sample" = "individual"))
  
atp1a1a.genotype <- ggplot(atp1a1a.df.geno) +
  geom_tile(aes(POS, sample, fill = factor(genotype, levels = c(2,1, 0)))) +
  # geom_point(aes(as.numeric(POS), sample, color = factor(genotype, levels = c(2,1,0)))) +
  scale_fill_discrete(palette =  c("#FFC107", "#D81B60","#1E88E5"), name =  "haplotype") +
  scale_color_discrete(palette =  c("#FFC107", "#D81B60","#1E88E5"), name =  "haplotype") +
  # scale_fill_gradient2(low = "#009E73", mid = "#E69F00", high = "#56B4E9", midpoint=0.5, name =  "MDS Scaled") +
  #scale_fill_gradient2(low = "firebrick3", mid = "orange" ,high = "darkgreen", midpoint=0.5, name =  "MDS Scaled") +
  #scale_x_continuous(labels = function(x) paste0(x / 1e6),name = "Mbps", expand = c(0.01,0)) +
  facet_grid(Ecotype+Population~.,scale = "free", space = "free", switch = "y") +
  theme_classic() +
  ggtitle("Atp1a1a t1") +
  theme(legend.position = "bottom", panel.spacing.y = unit(0,'lines'), panel.spacing.x = unit(0.5,'lines'),
        legend.frame = element_rect(colour="black"),
        legend.ticks = element_line(colour="black"),
        axis.title.y.right = element_blank(),                # hide right axis title
        axis.text.y.right = element_blank(),                 # hide right axis labels
        axis.ticks.y = element_blank(),                   # hide left/right axis ticks  
        # axis.text.y = element_blank(),                    # hide left/right axis ticks
        axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        # axis.text.y = element_text(margin = margin(r = 0), size = 5),  # move left axis labels closer to axis 
        strip.background = element_rect(color = "black", size = 0.5),
        panel.background = element_rect(fill = "grey90", color = "black", size = 0.5))


#  Calc nj 
library(ape)
library(ggtree)
library(patchwork)
nj.data <- nj(dc.atp1a1a)

# Create tree plot
plot.tree <- ggtree(nj.data, aes(color = Ecotype), layout = "ape")
#  Combine with sample data
plot.tree <- plot.tree %<+% mds.atp1a1a.df
## Custom tip colours
plot.tree <- plot.tree + geom_tippoint(aes(color = Ecotype, shape = Ecotype), stroke = 1, size=3) +
  geom_tiplab(aes(label = paste(Population, label)),hjust = -0.25, size = 2) +
  scale_shape(na.translate = TRUE, na.value = 21) + 
  scale_x_continuous(expand = c(0.25,0.25)) +
  scale_y_continuous(expand = c(0.25,0.25)) +
  ggtitle("Atp1a1a t1") +
  theme(legend.position = "bottom") 

ggsave(paste0(outdir,"/atp1a1a_t1_genotype.png"), 
      atp1a1a.genotype + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20, height = 12)
#ggsave(paste0("test.png"), atp1a1a.genotype + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20, height = 12)
#ggsave(paste0("test.pdf"), atp1a1a.genotype + plot.tree + (p.het.atp1a1a/p.mds.atp1a1a) + plot_layout(widths = c(1,1,1)), width = 20, height = 12)


## Calculate nj tree for all genes

gene.df <- as.data.table(gene.bed[gene.bed$type=="gene",])

gene.calcs <- map2_dfr(
    gene.df$chr, seq_len(nrow(gene.df)),
    ~{
      chr <- .x
      idx <- .y
      # chr <- "CM102076.1"
      # idx <- 20
      wnd <- gene.df[idx]
      wndname <- paste0(chr,"-", wnd$start,"-",wnd$end, "-", wnd$gene)


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
      if(ncol(geno.mat)<=20|nrow(geno.mat)<=20){ return(NULL) ; print("Not enough SNPs or samples")}
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
      nj.tree <- nj(dc)
      ## Get nj tree
      nj.tree <- nj(dc)

      mds.df <- data.frame(sample = colnames(geno.mat), chr = chr, start = wnd$start, end = wnd$end, Duke.gene = wnd$Duke.gene, gene = wnd$gene, base.calls = base.calls, 
                        nsamps = ncol(geno.mat), nsnps = nrow(geno.mat),
                        MDS1 = mds[,1], MDS2 = mds[,2], dist.marine = dc.marine.samp, het.prop = het.count/base.calls, tree = write.tree(nj.tree))

      

      # Return data.frame with PCA1, PCA2, MDS1, MDS2, sample names, window name
      return(mds.df)

    }
  )

unique(gene.calcs$gene)

# Merge with sample data
gene_data <- merge(gene.calcs, seq_data,  by.x = "sample", by.y="individual", all.x = T) %>%
                  select(-tree) %>%
                  mutate(windowname = paste(chr, start, end, sep = "-"),
                         mid = end - ((end-start)/2))

## Extract trees (remove duplicated trees)
gene.calcs.tree <- gene.calcs[!duplicated(gene.calcs$Duke.gene),] %>%
        arrange(start)

gene.mds.plot <- ggplot(gene_data) +
  geom_jitter(aes(MDS1, start, col = Ecotype)) +
  geom_text(aes(min(MDS1*1.1), start, label = gene))

ggsave("test.png", gene.mds.plot)


# Function that rotates trees so they are polarised
rotate_points <- function(x, y, angle_deg, center = c(0, 0)) {
  # Convert degrees to radians
  angle_rad <- pi * (angle_deg / 180)
  
  # Shift to center
  x_shift <- x - center[1]
  y_shift <- y - center[2]
  
  # Rotation matrix
  cos_a <- cos(angle_rad)
  sin_a <- sin(angle_rad)
  
  x_rot <- x_shift * cos_a - y_shift * sin_a + center[1]
  y_rot <- x_shift * sin_a + y_shift * cos_a + center[2]
  
  data.frame(x_rot = x_rot, y_rot = y_rot)
}

# range of angle to test rot 
angles_range <-seq(-180, 180, 5)

# Loop through and plot
nj.plots <- list()
for(i in 1:nrow(gene.calcs.tree)){
  # Name of window
  gene.tmp <- ifelse(is.na(gene.calcs.tree$gene[i]), gene.calcs.tree$Duke.gene[i], gene.calcs.tree$gene[i])
  # Window infor
  # plot
  plot.tmp <- ggtree(read.tree(text = gene.calcs.tree$tree[i]), layout = "ape") 
  # Combine with sample info
  plot.tmp <- plot.tmp %<+% seq_data
  
  # What is the middle coordinates of tree
  mid.point <- c((max(plot.tmp$data$x)+min(plot.tmp$data$x))/2, (max(plot.tmp$data$y)-min(plot.tmp$data$y))/2)
  
  # If you rotate the mig samples at what point is their x height min. (How do you rotate the tree to get the mig samples at the bottom)
  best_rot <- angles_range[which.min(sapply(angles_range, FUN = function(x){
  mean(rotate_points(plot.tmp$data$x[plot.tmp$data$Ecotype=="mig"], plot.tmp$data$y[plot.tmp$data$Ecotype=="mig"], x, center = mid.point)[,"y_rot"], na.rm = T)
  }))]
    
  plot.tmp$data[,c("x", "y")] <- rotate_points(plot.tmp$data$x+min(plot.tmp$data$x), plot.tmp$data$y, best_rot, center = mid.point)
  plot.tmp$data[,c("branch.x", "branch.y")] <- rotate_points(plot.tmp$data$branch.x, plot.tmp$data$branch.y, best_rot, center = mid.point)
  
  plot.tmp.rot <- ggplot(plot.tmp$data) +
    # geom_point(aes(x_new, y_new)) +
    # geom_point(aes(mid.point[1], mid.point[2]), col = "red") +
    geom_segment(aes(x = x, y=y, xend = branch.x+(branch.x-x), yend = branch.y+(branch.y-y))) +
    geom_point(aes(x, y, shape = Ecotype, col = Ecotype), size = 2) +
    scale_color_manual(values = c("#FFC107", "#1E88E5", "#009E73"), na.translate = F) +
    scale_shape_manual(values = c(15,19,17), na.translate = F) +
    # scale_color_manual(values = c("#E69F00", "#009E73","#D55E00","#0072B2")) +
    coord_fixed() +
    theme_void() +
    theme(text = element_text(size = 20))

    # Alternate title top and bottom
    if((i %% 2) == 1){plot.tmp.rot <- plot.tmp.rot + ggtitle(gene.tmp) + theme(plot.title = element_text(hjust = 0.5, size = 15))}
    if((i %% 2) == 0){plot.tmp.rot <- plot.tmp.rot + labs(caption = gene.tmp) + theme(plot.caption = element_text(hjust=0.5, size = 15))}

  # Add tree tips
  nj.plots[[i]] <- plot.tmp.rot
  print(paste("Done:", gene.tmp, i))
}

# 1. Standard creation of the combined sub-plots
nj.gene.plots.comb <- wrap_plots(nj.plots, guides = "collect", widths = 1, heights = 1, nrow = 1) & 
                      theme(legend.position = "right")

# 2. Treat the inner assembly as a single element using wrap_elements()
# This forces patchwork to give the entire row just ONE tag (b) 
nj.gene.plots.comb_wrapped <- wrap_elements(nj.gene.plots.comb)

# 3. Combine your top plot with the wrapped lower element and save
ggsave("test.png", 
       (mds_chrI_inv_genes / nj.gene.plots.comb_wrapped + plot_layout(heights = c(6, 2))) + 
       plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.tag = element_text(size = 16))), 
       width = 15, height = 16)
       
ggsave(paste0(outdir,"/chrI_inv_mds_and_gene-njtree.png"),
       (mds_chrI_inv_genes / nj.gene.plots.comb_wrapped + plot_layout(heights = c(6, 2))) + 
       plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.tag = element_text(size = 16))), 
       width = 15, height = 16)

ggsave(paste0(outdir,"/chrI_inv_mds_and_gene-njtree.tiff"),
       (mds_chrI_inv_genes / nj.gene.plots.comb_wrapped + plot_layout(heights = c(6, 2))) + 
       plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.tag = element_text(size = 16))), 
       width = 15, height = 16)

ggsave(paste0(outdir,"/chrI_inv_mds_and_gene-njtree.pdf"),
       (mds_chrI_inv_genes / nj.gene.plots.comb_wrapped + plot_layout(heights = c(6, 2))) + 
       plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.tag = element_text(size = 16))), 
       width = 15, height = 16)

### Fst ###
library(hierfstat)

# Subset just to Uist samples
Uist.samples <- seq_data$individual[seq_data$Region=="Uist"]
Uist.samples <- Uist.samples[Uist.samples%in%(colnames(vcf.atp1a1a@gt)[-1])]

# convert to genind
my_genind.atp1a1a <- vcfR2genind(vcf.atp1a1a[,,samples = which((colnames(vcf.atp1a1a@gt)[-1])%in%Uist.samples)], sep = "/", return.alleles = TRUE)
my_genind.atp1a1a@pop <- as.factor(seq_data$Ecotype[match(rownames(my_genind.atp1a1a@tab), seq_data$individual)])

# Calculate Fst
pop.diff <- genet.dist(my_genind.atp1a1a, diploid = T, method = "WC84")
pop.diff
# Save results
write.table(file = paste0(outdir,"/atp1a1a_t1_Fst.txt"), x = as.matrix(pop.diff))

### Write out phyfile
#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes,
# Extract atp1a1a genotype
gt.apt1a1a <- extract.gt(vcf.atp1a1a[is.polymorphic(vcf.atp1a1a, na.omit = T),], element = "GT", return.alleles = T)
# Remove non Uist samples
gt.apt1a1a <- gt.apt1a1a[,colnames(gt.apt1a1a)%in%(seq_data$individual[seq_data$Region=="Uist"])]
# Remove high heterozygous individuals
gt.apt1a1a <- gt.apt1a1a[,!colnames(gt.apt1a1a)%in%high.het.ind]
# Remove high missing samples
gt.apt1a1a <- gt.apt1a1a[,!colnames(gt.apt1a1a)%in%high.cov.ind]


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
ape::write.dna(t(gt.apt1a1a), file = paste0(outdir,"/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.homoInd.phy"), format ="interleaved")