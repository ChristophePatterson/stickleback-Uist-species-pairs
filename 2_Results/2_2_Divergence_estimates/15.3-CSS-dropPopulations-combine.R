## devtools::install_github("jdstorey/qvalue")
library(qvalue)
library(tidyverse)

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

args<-commandArgs(trailingOnly=T)
CSS.dir <- args[1]
CSS.dir <-  "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/sliding-window/CSS/dropPops"

#CSS.run <- "stickleback.wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05.CSSm.10000perm.txt"

## Get combinations of dropped pops
dropComb <- data.frame(waterbody = c("CLAC", "DUIN", "LUIB", "OBSE"))

dropComb$dropped <- paste0("n", dropComb$waterbody)
# Getfile locations
dropComb$run_name <- paste0("stickleback.", dropComb$dropped, ".wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05")
dropComb$run_dir <- paste0(CSS.dir, "/", dropComb$dropped, "/", dropComb$run_name, "/")
dropComb$run_file <- paste0(dropComb$run_dir, dropComb$run_name, ".CSSm.10000perm.txt")
dropComb$run_file

## Read in chrom info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))

## Read in each CSS calc
CSS.list <- apply(dropComb, 1, function(x) {
    df <- cbind.data.frame(
      dropped = x["dropped"],
      read.table(x["run_file"], header = TRUE)
    )
    df <- df[df$nsnps>=5|df$nsnps<=200,]
    df$qval.0001 <- qvalue(1 - df$pval, fdr.level = 0.0001)$qvalues
    df
  })

# Combined into single data frame
CSS.long <- do.call("rbind", CSS.list)

CSS.long$chr <- chr$Sequence.name[match(CSS.long$chr, chr$GenBank.seq.accession)]
CSS.long$chr <- factor(CSS.long$chr, levels = chr$Sequence.name)

CSS.long <- CSS.long %>%
    mutate(qval.sig.0001 = qval.0001 < 0.05)

CSS.long <- CSS.long %>%
  group_by(dropped) %>%
  mutate(temp = replace(qval.sig.0001, is.na(qval.sig.0001), FALSE), 
         temp1 = cumsum(!temp)) %>%
  group_by(temp1, .add = T) %>%
  mutate(goal =  +(row_number() == which.max(temp) & any(temp))) %>%
  group_by(chr) %>%
  group_by(dropped, .add = T) %>%
  mutate(goal.0001 = ifelse(temp, paste0(first(dropped), "_", first(chr), "_", cumsum(goal)), NA)) %>%
  select(-temp, -temp1)

unique(CSS.long$goal.0001)
top.regions <- sort(table(CSS.long$goal.0001[CSS.long$css>=2]), decreasing = T)*500
top.regions[top.regions>=2500]

# Pivot table wider
CSS.wide <- CSS.long %>%
  pivot_wider(
    names_from = dropped,
    values_from = c(nsnps, css, nperms, pval, qval.sig.0001, qval.0001, goal.0001)
  )

png("test.png", width = 1000, height = 1000)
plot(CSS.wide[,grepl("css", colnames(CSS.wide))])
dev.off()

p <- ggplot(CSS.long) +
    geom_point(aes(start, css, col = qval.sig.0001)) +
    geom_point(data = CSS.long[CSS.long$goal.0001%in%names(top.regions[top.regions>=2500]),], aes(start, css), col = "red") +
    facet_grid(dropped~chr, scale = "free_x")

ggsave("test.png", p, width = 20, height = 10)

p <- ggplot(CSS.long) +
    geom_segment(aes(x = start, xend = end, y = dropped, yend= dropped, col = qval.sig.0001), linewidth = 3) +
    geom_segment(data = CSS.long[CSS.long$goal.0001%in%names(top.regions[top.regions>=2500]),], aes(x = start, xend = end, y = dropped, yend= dropped), col = "red", linewidth = 3) +
    facet_grid(chr~., scale = "free_x")

ggsave("test.png", p, width = 8, height = 15)

