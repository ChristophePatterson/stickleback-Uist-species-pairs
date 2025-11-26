## Chromosome 1 sliding window plot

#libraries
library(vcfR)
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(patchwork)

# Christmas  color
xmas = c(
  "gold", # gold
  "#78DFE4", # Sky Blue Crayola
  "#C8001C", # Lava
  "#7C0322", # Burgundy
  "#91C153", # Pistachio
  "#609FA4"  # Cadet Blue
)

##### ## Read in vcf
##### vcf.SNPs <- read.vcfR("/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX.vcf.gz")
##### 
##### ## Waterbodies of interest
##### paired_sp_waterbodies <- c("DUIN", "OBSE", "LUIB", "CLAC", "OLAV", "TORM")
##### 
##### ### Combine Uist samples with Uist sample locations
##### seq_data <- as_tibble(read.csv("bigdata_Christophe_header_2025-04-28.csv", header = T))
##### 
##### seq_data$Ecotype[seq_data$Ecotype=="st"] <- "fw"
##### seq_data$Ecotype[seq_data$Ecotype=="anad"] <- "mig"
##### 
##### # seq_data <- seq_data[seq_data$Ecotype!="fw",]
##### 
##### ## REmove unneeded samples
##### seq_data <- seq_data[seq_data$Waterbody%in%paired_sp_waterbodies,]
##### samples <- data.frame(samples = colnames(vcf.SNPs@gt)[-1])
##### unique(seq_data$Waterbody)
##### vcf.SNPs <- vcf.SNPs[samples = samples$samples[samples$samples%in%seq_data$individual]]
##### 
##### # Remove invaient and non-biallelic snps
##### vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
##### vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs, na.omit = T),]
##### 
## Get list of chromosomes
chr <- unique(vcf.SNPs@fix[,"CHROM"])

# Create windows 
wndsize <- 1000000
wndslid <- 1000000

# Create data set of window positions
chr.df <- as.data.table(vcf.SNPs@fix)[, .(POS = as.numeric(POS)), by = CHROM]
sldwindows.df <- chr.df[, {
  end <- seq(min(POS), max(POS), by = wndslid)
  start <- end - wndsize + 1
  .(start = start, end = end)
}, by = CHROM]


# Calculate number of snps within each window
sldwindows.df$nsnps <- apply(sldwindows.df, MARGIN = 1, function(x) sum(chr.df$CHROM==x["CHROM"] & between(as.numeric(chr.df$POS), x["start"], x["end"])))

# Add in sliding window chr name
# Read in chr name df
scaf <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
sldwindows.df$chr_name <- scaf$Sequence.name[match(sldwindows.df$CHROM, scaf$GenBank.seq.accession)]

# Extract geno
geno.mat.full <- extract.gt(vcf.SNPs, element = "GT")

# Remove windows with less than 5 snps
sldwindows.df <- sldwindows.df[sldwindows.df$nsnps>=100,]

print("Starting PCA and MDS run")

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

## Loop through each window
nj.sldwnd <- list()
nj.sldwnd <- apply(sldwindows.df, MARGIN = 1, function(x) {
  # Extract window info
  wnd <- x
  # Get name for window
  wndname <- paste0(wnd["chr_name"],"-", wnd["start"],"-",wnd["end"], "-n", wnd["nsnps"])
  
  # Get snps ID that are within window
  snp.idx <- which(vcf.SNPs@fix[, "CHROM"] == wnd["CHROM"] &
                     between(as.numeric(vcf.SNPs@fix[, "POS"]), wnd["start"], wnd["end"]))
  
  print(paste(wndname, "nsnps",length(snp.idx)))
  # If less than two snps remove
  if(length(snp.idx) < 2) { return(NULL) ; print("Not enough SNPs or samples")}
  ## Set geno to be 0,1, or 9
  # Extract Snps within window
  geno.mat <- geno.mat.full[snp.idx, , drop = FALSE]
  # Convert to matrix of 0,1, and 2
  geno.mat <- matrix(as.integer(factor(geno.mat, levels = c("0/0", "0/1", "1/1")))-1, 
                     nrow = nrow(geno.mat), dimnames = dimnames(geno.mat))
  geno.mat <- apply(geno.mat, 2, as.integer)
  
  # Filter rows/cols with excessive NA, constant values, etc.
  is.samp.all.missing <- apply(geno.mat, MARGIN = 2, function(x) sum(is.na(x))>ncol(geno.mat)*0.5)
  if(any(is.samp.all.missing)){geno.mat <- geno.mat[,-is.samp.all.missing]}
  
  is.bad.snp <- apply(geno.mat, 1, function(x) {
    vals <- unique(na.omit(x))
    length(vals) <= 1 || all(vals == 1)  # All het or constant
  })
  if(any(is.bad.snp)){geno.mat <- geno.mat[-which(is.bad.snp), ]}
  # print(paste("Stage 1", nrow(geno.mat), ncol(geno.mat)))
  #Check is matrix still contains data
  if(ncol(geno.mat)<=10|nrow(geno.mat)<=10){ return(list(wndname, NA)) ; print("Not enough SNPs or samples")}
  # Run MDS / PCA on cleaned geno.sub
  # print("Stage 2")
  ## Calculate distance matrix
  dc <- dist(t(geno.mat))
  # Check is dist matrix contains any NA
  if(any(is.na(dc))) {return(list(wndname, NA)) ; print("mds-failed")}
  # Calculate nj tree
  nj.tree <- nj(dc)
  
  # Return 
  return(list(wndname, nj.tree))
  
}
)

# Remove list elements that are NA
nj.sldwnd <- nj.sldwnd[!is.na(sapply(nj.sldwnd, "[[", 2))]
# Remove list elements that are NULL
nj.sldwnd <- nj.sldwnd[!sapply(nj.sldwnd, is.null)]
# Location of atp1a1
ATP1A1 <- data.frame(chr = "I", start = 26848100, end = 26861753, name = "atp1a1a", Ecotype="", Population="")
i <- 1

# range of angle to test rot 
angles_range <-seq(-180, 180, 5)

# Loop through and plot
nj.plots <- list()
tree.data.list <- list()
for(i in 1:length(nj.sldwnd)){
  # Name of window
  wndname.tmp <- nj.sldwnd[[i]][[1]]
  # Window infor
  wnd.prop.tmp <- c(stringr::str_split(wndname.tmp, "-"))
  # plot
  plot.tmp <- ggtree(nj.sldwnd[[i]][[2]], layout = "ape") 
  # Combine with sample info
  plot.tmp <- plot.tmp %<+% seq_data
  
  mid.point <- c((max(plot.tmp$data$x)+min(plot.tmp$data$x))/2, (max(plot.tmp$data$y)-min(plot.tmp$data$y))/2)
  
  # Find best rotation to have migrants at the bottom
  best_rot <- angles_range[which.min(sapply(angles_range, FUN = function(x){
    median(rotate_points(plot.tmp$data$x[plot.tmp$data$Ecotype=="mig"], plot.tmp$data$y[plot.tmp$data$Ecotype=="mig"], x, center = mid.point)[,"y_rot"], na.rm = T)
  }))]
  
  
  plot.tmp$data[,c("x", "y")] <- rotate_points(plot.tmp$data$x+min(plot.tmp$data$x), plot.tmp$data$y, best_rot, center = mid.point)
  plot.tmp$data[,c("branch.x", "branch.y")] <- rotate_points(plot.tmp$data$branch.x, plot.tmp$data$branch.y, best_rot, center = mid.point)
  
  
  plot.tmp.rot <- ggplot(plot.tmp$data) +
    # geom_point(aes(x_new, y_new)) +
    # geom_point(aes(mid.point[1], mid.point[2]), col = "red") +
    geom_segment(aes(x = x, y=y, xend = branch.x+(branch.x-x), yend = branch.y+(branch.y-y))) +
    geom_point(aes(x, y, shape = Ecotype, col = Ecotype)) +
    scale_color_manual(values = xmas, na.value = "NA") +
    coord_fixed() +
    theme_void()
  
  # Add tree tips
  nj.plots[[i]] <- plot.tmp.rot + ggtitle(gsub(" ","",nj.sldwnd[[i]][[1]]))
  print(paste0("Done: ", nj.sldwnd[[i]][[1]], i))

  # Extract data
  plot.data.tmp <- plot.tmp$data[,c("node", "parent", "x", "y","branch.x", "branch.y", "isTip", "label", "Ecotype", "Population", "Waterbody")]
  plot.data.tmp$window <- wndname.tmp
  plot.data.tmp$CHROM <- wnd.prop.tmp[[1]][1]
  plot.data.tmp$start <- as.numeric(wnd.prop.tmp[[1]][2])
  plot.data.tmp$end <- as.numeric(wnd.prop.tmp[[1]][3])
  tree.data.list[[i]] <- plot.data.tmp

}

# Exmple plots
# ggsave("test.png", plot = wrap_plots(nj.plots[[1]], guides = "collect", widths = 1, heights = 1) , width = 20, height = 20)

# Merge plots
nj.plots.comb <- wrap_plots(nj.plots, guides = "collect", widths = 1, heights = 1) & theme(legend.position = "bottom", title = element_text(size = 5))

# Save test
ggsave("test.png", plot = nj.plots.comb, width = 40, height = 40)
# Save
ggsave(paste0("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/njtree/njtree_wnd",wndslid,"sld",wndslid,".png"), plot =nj.plots.comb, width = 40, height = 40)

# Combine tree data
tree.data.df <- do.call(rbind, tree.data.list)

## Change order of chromosome based on position
tree.data.df$CHROM <- factor(tree.data.df$CHROM, levels = scaf$Sequence.name)

p <- ggplot(tree.data.df) +
  geom_segment(aes(x = x, y=y, xend = branch.x+(branch.x-x), yend = branch.y+(branch.y-y)), linewidth = 0.2) +
  geom_point(aes(x, y, shape = Ecotype, col = Ecotype), size =  0.5) +
  scale_color_manual(values = xmas, na.value = "NA") +
  coord_fixed() +
  theme_void() + 
  facet_grid(CHROM~round(start/wndsize, 1), switch = "y") +
  # Make background white
  theme(plot.background = element_rect(fill = "white", color = "white"),
  # Make y facet labels visible
        strip.text.y = element_text(face = "bold", angle = 90),
                text = element_text(size = 25)
  ) 

# ggsave("test.png", plot = p, width = 40, height = 30)
# Save
ggsave(paste0("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/njtree/njtree_wnd",wndslid/1000,"Kb_sld",wndslid/1000,"Kb_bychr_scaled.png"), plot = p, width = 40, height = 30)
