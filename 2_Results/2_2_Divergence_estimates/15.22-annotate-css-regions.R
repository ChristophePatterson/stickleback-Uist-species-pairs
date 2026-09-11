library(tidyverse)
library(patchwork)
library(IRanges)
library(ggrepel)

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

## args<-commandArgs(trailingOnly=T)
## CSS.dir <- args[1]
CSS.dir <-  "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops"
Genome.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/genomes"

CSS.run <- ".wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05"

top.grouped.regions.table <-  read_csv(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_all_sig_top_regions_grouped.10Kbp.txt"))

top.grouped.regions.table %>%
    arrange(desc(mn.CSS.sig))

## Read in chromsome file

## Read in chrom info
genome.info <- as_tibble(read.table(paste0(Genome.dir, "/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv"), sep = "\t", header = T))
genome.info$Sequence.name <- gsub("chr", "", genome.info$Sequence.name)
# Calculate cummulative sum
genome.info$Cum.Seq.length <- c(0, cumsum(genome.info$Seq.length[1:(nrow(genome.info)-1)]))
chr.order <- genome.info$Sequence.name

# Read in CSS droppop results
CSS.long <- read_csv(paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_combine.csv")) %>%
  mutate(chr = factor(chr, levels = chr.order))

# Read in DUKE annotation
DUKE.gtf <- read_table(paste0(Genome.dir, "/GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf"), col_names = T)

DUKE.gtf$chr <- genome.info$Sequence.name[match(DUKE.gtf$chr, genome.info$GenBank.seq.accession)]

# Pivot table wider
CSS.wide <- CSS.long %>%
  pivot_wider(
    names_from = dropped,
    values_from = c(nsnps, css, nperms, pval, qval.sig.0001, qval.0001, goal.0001)
  )

## Calculate cumulative start and end point
# Start
CSS.wide$start.cum <- CSS.wide$start+(genome.info$Cum.Seq.length[match(CSS.wide$chr,genome.info$Sequence.name)])
DUKE.gtf$start.cum <- DUKE.gtf$start+(genome.info$Cum.Seq.length[match(DUKE.gtf$chr,genome.info$Sequence.name)])
# v5.mpDUKE.gtf$start.cum <- v5.mpDUKE.gtf$start+(genome.info$Cum.Seq.length[match(v5.mpDUKE.gtf$chr,genome.info$Sequence.name)])
# End
CSS.wide$end.cum <- CSS.wide$end+(genome.info$Cum.Seq.length[match(CSS.wide$chr,genome.info$Sequence.name)])
DUKE.gtf$end.cum <- DUKE.gtf$end+(genome.info$Cum.Seq.length[match(DUKE.gtf$chr,genome.info$Sequence.name)])
# v5.mpDUKE.gtf$end.cum <- v5.mpDUKE.gtf$end+(genome.info$Cum.Seq.length[match(v5.mpDUKE.gtf$chr,genome.info$Sequence.name)])

# Calculate mean across all populations
CSS.wide$mn.CSS <- apply(CSS.wide, MARGIN = 1, function(x) mean(as.numeric(x[grep("css_", names(x))]), na.rm = T))
CSS.wide$mx.CSS <- apply(CSS.wide, MARGIN = 1, function(x) max(as.numeric(x[grep("css_", names(x))]), na.rm = T))

## Determine which windows are significant across all dropped population comparisons
CSS.wide$all.sig.qvalue.0001 <- apply(CSS.wide[,grepl("qval.sig.0001_", colnames(CSS.wide))], MARGIN = 1, function(x) 
  ifelse(any(is.na(x)), return(NA), all(x)))

## Determine which windows are significant across any dropped population comparisons
CSS.wide$any.sig.qvalue.0001 <- apply(CSS.wide[,grepl("qval.sig.0001_", colnames(CSS.wide))], MARGIN = 1, function(x) any(x))

## Find all the contigous regions
CSS.wide <- CSS.wide %>%
  mutate(temp = replace(all.sig.qvalue.0001, is.na(all.sig.qvalue.0001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1, .add = T) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  mutate(all.goal.0001 = ifelse(temp, paste0(dplyr::first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

## Table
CSS.wide$any.sig.qvalue.0001
top.regions.table <- CSS.wide %>%
  filter(!is.na(all.goal.0001)) %>%
  group_by(all.goal.0001) %>%
  summarise(chr = dplyr::first(chr),
            start = min(start),
            end = max(end),
            mn.CSS.sig = mean(mn.CSS[all.sig.qvalue.0001&(!is.na(all.sig.qvalue.0001))]), # Must be above the mn.CSS calc
            mn.CSS = mean(mn.CSS),
            mx.CSS = max(mx.CSS),
            wnd.length = (1+max(end))-min(start)) %>%
  ungroup() %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(bp.break = c(Inf, diff(start)))


## Interesting regions
regions <- data.frame(chr = factor(c("I", "IV", "IX", "XI", "XXI"), levels = genome.info$Sequence.name), start = c(26000000, 12000000, 12000000, 5500000, 9000000), end = c(27500000, 16000000, 14500000,7000000, 12500000))

## top regions
top.regions.table %>%
  arrange(desc(mn.CSS)) %>%
  print(n = 50)

# Subset to top 10 regions 
regions <- top.grouped.regions.table %>%
    arrange(desc(mn.CSS.sig)) %>%
    filter(row_number() <= 10) %>%
    select(chr, start, end) 

# Custom regions to plot
regions <- data.frame(chr = factor(c("I",     "IV",     "IV",      "IV",    "XIV",   "XIX"), levels = genome.info$Sequence.name), 
                      start = c(26500000, 12820000, 27500000,  27500000, 12900000, 2740000), 
                      end   = c(27200000, 13050000, 27900000,  27800000, 13000000, 2880000))

regions$start.cum <- regions$start+(genome.info$Cum.Seq.length[match(regions$chr, genome.info$Sequence.name)]+genome.info$Cum.Seq.length[1])

genome.all.plot <- ggplot(genome.info) +
  geom_segment(aes(x = Cum.Seq.length, xend = Cum.Seq.length+Seq.length, y = 0, col = as.numeric(Chromosome.name) %% 2 == 0), linewidth = 10) +
  geom_segment(data = regions, aes(x = start.cum, xend = start.cum+(start-end), y = 0), linewidth = 10, col = "firebrick") +
  geom_point(data = regions, aes(x = start.cum+(start-end)/2, y = 0), shape = 4, size = 3, col = "firebrick", stroke = 1) +
  geom_text(aes(Cum.Seq.length+(Seq.length/2), y = 0, label = Sequence.name), col = "white") + 
  # geom_vline(aes(xintercept = Cum.Seq.length)) +
  scale_color_manual(values = c("grey30", "grey50")) +
  scale_y_continuous(expand = c(0.05,0.05)) +
  scale_x_continuous(expand = c(0,0), label = function(x) x/1e6, name = "Mbps") +
  theme_void() +
  theme(legend.position = "none")
genome.all.plot

# Create a list to store the chromosome plots
plot.chr <- list()

i <- 1
for(i in 1:nrow(regions)){
    # Get chromosome, start and end for the region
  chr.i <- regions$chr[i]
  length.i <- regions$end[i] - regions$start[i]
  ## Optional if length is less than 100Kb, then set length to 100Kbp
  # if(length.i < 100000){length.i <- 100000}
  start.i <- regions$start[i]  # - (length.i*0.5)
  end.i <- regions$end[i]      # + (length.i*0.5)
  # If start is less than 0, set start to 0
  if(start.i < 0){start.i <- 0}

  # Loop print statement to check which region is being plotted
  print(paste0("Plot window for ", chr.i, ": ", start.i,"-",end.i))
  # Subset CSS data to the region
  CSS.temp <- CSS.long[CSS.long$chr==chr.i&between(CSS.long$start, start.i, end.i),]

  # Create plot for the region
  plot.chr[[i]] <- ggplot(CSS.temp) +
  geom_line(aes(start, css, col = dropped, group = dropped), linewidth = 1) +
  geom_segment(data = DUKE.gtf[DUKE.gtf$chr==chr.i&DUKE.gtf$type=="gene",], 
               aes(x = start, xend = end, y = -0.4), linewidth = 3, col = "cornflowerblue") +
  geom_segment(data = DUKE.gtf[DUKE.gtf$chr==chr.i&DUKE.gtf$type=="exon",], 
               aes(x = start, xend = end, y = -0.4), linewidth = 3, col = "orange") +
  geom_text(data = DUKE.gtf[DUKE.gtf$chr==chr.i&DUKE.gtf$type=="gene",], 
              aes(x = end-((end-start)/2), y = -0.8, label = gene_id), size = 2) +
  geom_text(data = DUKE.gtf[DUKE.gtf$chr==chr.i&DUKE.gtf$type=="gene",], 
               aes(x = end-((end-start)/2), y = -1.2, label = fGas.name.match), size = 2) +
  geom_text(data = DUKE.gtf[DUKE.gtf$chr==chr.i&DUKE.gtf$type=="gene",], 
            aes(x = end-((end-start)/2), y = -1.6, label = v5.name.match), size = 2) +
  scale_color_manual(values = cbPalette) +
  # limit x axis to the region but still draw segement that start outside of region
  # don't remove segments that overlap with plot boundary ggplot
  coord_cartesian(xlim = c(start.i, end.i), expand = T) +
  ggtitle(paste0(chr.i, ": ", start.i,"-",end.i)) +
  labs(x = "Genomic Position (bp)", y = "CSS") + 
  theme_bw()
  
}

# Combine all chromosome plots with the genome plot
plot.chr.all <- (genome.all.plot / wrap_plots(plot.chr, ncol = 1,)) + plot_layout(heights = c(1, 20))
# Save plot
ggsave(paste0("test.pdf"), plot.chr.all, height = 30, width = 20)
ggsave(paste0(CSS.dir, "/CSS_zoom_annotated.pdf"), plot.chr.all, height = 20, width = 20)
ggsave(paste0(CSS.dir, "/CSS_zoom_annotated.png"), plot.chr.all, height = 20, width = 20)

# Convert to IGV format
# subset CSS to regions of interest
CSS.wide.filt <- CSS.wide %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

# Create IGV format table
CSS.wide.igv <- CSS.wide.filt %>%
  select(chr, start, end, mn.CSS) %>%
  mutate(start = as.integer(start), end = as.integer(end), mn.CSS = as.numeric(mn.CSS)) %>%
  arrange(chr, start)
# Overwrite chromsome name with GenBank accession for IGV
CSS.wide.igv$chr <- genome.info$GenBank.seq.accession[match(CSS.wide.igv$chr, genome.info$Sequence.name)]
# Write table
write.table(CSS.wide.igv, paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_CSS_combine_igv.BEDgraph"), sep = "\t", quote = F, row.names = F, col.names = F)  

# Subset DUKE annotation to regions of interest
DUKE.gtf.filt <- DUKE.gtf %>%
  rowwise() %>%
  filter(any(
    chr == regions$chr &
      start >= regions$start &
      end <= regions$end
  )) %>%
  ungroup()

# Create IGV format table
DUKE.gtf.igv <- DUKE.gtf.filt %>%
  select(chr, start, end, fGas.name.match, v5.name.match, gene_id, method, type) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  arrange(chr, start)
# Overwrite chromsome name with GenBank accession for IGV
DUKE.gtf.igv$chr <- genome.info$GenBank.seq.accession[match(DUKE.gtf.igv$chr, genome.info$Sequence.name)]
# Write table
write.table(DUKE.gtf.igv, paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_DUKE_annotation_igv.gft"), sep = "\t", quote = F, row.names = F, col.names = F)

## Save regions of interest as BED file for IGV
regions.igv <- regions %>%
  select(chr, start, end) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  arrange(chr, start)
# Overwrite chromsome name with GenBank accession for IGV
regions.igv$chr <- genome.info$GenBank.seq.accession[match(regions.igv$chr, genome.info$Sequence.name)]
# Write table
write.table(regions.igv, paste0(CSS.dir, "/stickleback.dropPops.", CSS.run,"_regions_of_interest_igv.BED"), sep = "\t", quote = F, row.names = F, col.names = F)

# Read in gtf for duke
DUKE.gft.raw <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1/Duke_GAcu_1.gtf", header = F, sep = "\t")
DUKE.gft.raw$V1 <- genome.info$GenBank.seq.accession[match(as.numeric(gsub("chr", "", DUKE.gft.raw$V1)), genome.info$Chromosome.name)]
## Write table
write.table(DUKE.gft.raw, paste0("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1/Duke_GAcu_1_GenBankID.gtf"), sep = "\t", quote = F, row.names = F, col.names = F)
