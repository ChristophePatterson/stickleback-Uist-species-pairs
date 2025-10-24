## Load libraries
library(patchwork)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(data.table)
library(ggnewscale)

## Get arguments
args <- commandArgs(trailingOnly = TRUE)
fst_file_indPops <- args[1]
fst_file <- args[2]

## Read in Fst data for sliding windows between ecotypes
fst_file_indPops_df <- read_csv(paste0(fst_file_indPops, ".csv"))
## Read in Fst data for sliding windows between populations
fst_file_df <- read_csv(paste0(fst_file, ".csv"))

## Merge dataframes
GenPop_df <- fst_file_df %>%
  left_join(fst_file_indPops_df, by = c("scaffold","start","end","mid"), suffix = c("_pops","_ecotypes"))

## Label scaffolds as chromosomes
## Read in chromosome info
# Read in chr info
chr <- as_tibble(read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_sequence_report.tsv", sep = "\t", header = T))
## Calculate cumulative position
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))

# Replace scaffold name with chr name
GenPop_df$chr <- gsub("chr","",chr$Sequence.name[match(GenPop_df$scaffold, chr$GenBank.seq.accession)])
GenPop_df$chr <- factor(GenPop_df$chr, levels = gsub("chr", "", chr$Sequence.name[order(chr$GenBank.seq.accession)]))


# Calculate cumulative position of windows
GenPop_df$start.cum <- GenPop_df$start+(chr$Cum.Seq.length[match(GenPop_df$chr, gsub("chr", "", chr$Sequence.name))]+chr$Cum.Seq.length[1])

## Calculate which Fst are in the top 5% of the distribution
fst_threshold_95 <- quantile(GenPop_df$Fst_anad_resi, probs = 0.95, na.rm = T)
fst_threshold_99 <- quantile(GenPop_df$Fst_anad_resi, probs = 0.99, na.rm = T)
GenPop_df$anad_resi_outlier_95 <- ifelse(GenPop_df$Fst_anad_resi >= fst_threshold_95, "outlier", "non-outlier")
GenPop_df$anad_resi_outlier_99 <- ifelse(GenPop_df$Fst_anad_resi >= fst_threshold_99, "outlier", "non-outlier")

# Plot histogram of Fst values between ecotypes

p <- ggplot(GenPop_df[GenPop_df$Fst_anad_resi>=0,], aes(x=Fst_anad_resi)) +
  geom_histogram(binwidth = 0.01, fill="grey60", color="black") +
  geom_histogram(data=GenPop_df[GenPop_df$anad_resi_outlier_95=="outlier" & GenPop_df$Fst_anad_resi>=0,], 
        aes(x=Fst_anad_resi), binwidth = 0.01, fill="deepskyblue", color="black") +
  geom_histogram(data=GenPop_df[GenPop_df$anad_resi_outlier_99=="outlier" & GenPop_df$Fst_anad_resi>=0,], 
                 aes(x=Fst_anad_resi), binwidth = 0.01, fill="red", color="black") +
  scale_x_continuous(limits=c(0,0.25), expand=c(0,0)) +
  labs(x="Fst between ecotypes", y="Number of windows") +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", color = "black")
  )

q <- ggplot(GenPop_df[GenPop_df$Fst_anad_resi>=0,], aes(x=Fst_anad_resi)) +
  geom_histogram(binwidth = 0.01, fill="grey60", color="black") +
  geom_histogram(data=GenPop_df[GenPop_df$anad_resi_outlier_95=="outlier" & GenPop_df$Fst_anad_resi>=0,], 
                 aes(x=Fst_anad_resi), binwidth = 0.01, fill="deepskyblue", color="black") +
  geom_histogram(data=GenPop_df[GenPop_df$anad_resi_outlier_99=="outlier" & GenPop_df$Fst_anad_resi>=0,], 
                 aes(x=Fst_anad_resi), binwidth = 0.01, fill="red", color="black") +
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_bw() +
  labs(x="Fst between ecotypes", y="Number of windows") +
  coord_cartesian(ylim=c(0,200)) 

hist_plot <- q + inset_element(p,
  left = 0.25, bottom = 0.4, right = 0.95, top = 0.95,
  align_to = "panel",
  on_top = TRUE,
  clip = TRUE,
  ignore_tag = TRUE
)

ggsave(hist_plot, file="test.png", width=6, height=6)

### Create maxtrix of adaptive divergence Fst values between populations for outlier windows only
## Create long format dataframe
GenPop_df_long <- pivot_longer(GenPop_df[, c("chr","start","end","mid", "sites_pops", "sites_ecotypes", "anad_resi_outlier_95", "anad_resi_outlier_99" ,
                                colnames(GenPop_df)[grep("Fst_", colnames(GenPop_df))])], 
                                cols = colnames(GenPop_df)[grep("Fst_", colnames(GenPop_df))], values_to = "Fst", names_prefix = "Fst_", names_to = "Population")

## Average Fst per population pair for outlier windows only
adaptive_divergence_Fst <- GenPop_df_long %>%
    filter(Fst >=0, Population != c("anad_resi")) %>%
    group_by(Population) %>%
    group_by(anad_resi_outlier_95, .add = T) %>%
    summarise(mean_Fst = mean(Fst, na.rm = T),
            median_Fst = median(Fst, na.rm = T),
            n_windows = n(),
            total_sites_pops = sum(sites_pops, na.rm = T),
            total_sites_ecotypes = sum(sites_ecotypes, na.rm = T)) %>%
    mutate(Pop1 = str_split_i(Population, "_", 1),
            Pop2 = str_split_i(Population, "_", 2))

## Ensure Pop1 and Pop2 are factors with levels in correct order
## Order pop1 and pop2 so they are the same order as that in allpops
## Set order pops should be plotted
allpops <- c("CLAC" ,"DUIN","LUIB","OBSE", "CLAM", "DUIM", "LUIM", "OBSM")

## Order pop1 and pop2 so they are the same order as that in allpops
adaptive_divergence_Fst <- apply(adaptive_divergence_Fst, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops))]
}
      ) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(adaptive_divergence_Fst)

## Convert Pop1 and Pop2 to factors with correct levels
adaptive_divergence_Fst <- adaptive_divergence_Fst %>%
  mutate(Pop1 = factor(V1, levels = allpops),
         Pop2 = factor(V2, levels = allpops))

## Triangle for outlier windows only
outlier_triangle <- data.frame(x = c(0.4, 8.495, 0.4, 0.4), y = c(0.505, 8.6, 8.6, 0.6))
## Plot heatmap of mean Fst between populations for outlier windows only
Fst_outlier_comp_plot <- ggplot() +
    geom_tile(data = adaptive_divergence_Fst[adaptive_divergence_Fst$anad_resi_outlier_95=="non-outlier",], aes(x = Pop1, y = Pop2, fill = mean_Fst)) +
    geom_tile(data = adaptive_divergence_Fst[adaptive_divergence_Fst$anad_resi_outlier_95=="outlier",], aes(x = Pop2, y = Pop1, fill = mean_Fst)) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90", name = "Fst", limits = c(0, 0.3)) +
    theme_classic() +
    geom_vline(xintercept = 4.5) +
    geom_hline(yintercept = 4.5) +
    geom_hline(yintercept = 8.5) +
    geom_vline(xintercept = 8.5) +
    geom_polygon(data = outlier_triangle, aes(x, y), fill = "NA", color = "grey30", linewidth = 1.5) +
    geom_polygon(data = outlier_triangle, aes(y, x), fill = "NA", color = "deepskyblue", linewidth = 1.5) +
    theme(panel.background = element_rect(fill = "grey"), legend.frame = element_rect(colour = 'black'),
        axis.title = element_blank(), axis.text.y = element_text(angle=90, vjust = 4, hjust = 0.5),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        plot.margin = margin(b = 5.5, r = 5.5, t = 5.5, l = 30),
        axis.text.x = element_text(vjust=-1)) +
    annotate("text", x = -0.25, y = 6.5, label = "Migratory", angle = 90) + ## Axis annotation
    annotate("text", x = -0.25, y = 2.5, label = "Resident", angle = 90) +
    annotate("text", y = -0.25, x = 6.5, label = "Migratory") +
    annotate("text", y = -0.25, x = 2.5, label = "Resident") +
    annotate("text", y = 8.8, x = 2, label = "Non-outliers (95%)", col = "grey30") +
    annotate("text", y = 2, x = 8.8, label = "Outliers (5%)", angle = -90, col = "deepskyblue") +
    coord_fixed(clip = 'off', x = c(0.5, 8.5), y = c(0.5, 8.5))

## Save plot
ggsave(hist_plot + Fst_outlier_comp_plot, file=paste0(fst_file, "_adaptive_divergence_Fst_comparison_95perc.png"), width=12, height=6)

## Average Fst per population pair for outlier windows only
adaptive_divergence_Fst <- GenPop_df_long %>%
    filter(Fst >=0, Population != c("anad_resi")) %>%
    group_by(Population) %>%
    group_by(anad_resi_outlier_99, .add = T) %>%
    summarise(mean_Fst = mean(Fst, na.rm = T),
            median_Fst = median(Fst, na.rm = T),
            n_windows = n(),
            total_sites_pops = sum(sites_pops, na.rm = T),
            total_sites_ecotypes = sum(sites_ecotypes, na.rm = T)) %>%
    mutate(Pop1 = str_split_i(Population, "_", 1),
            Pop2 = str_split_i(Population, "_", 2))

## Ensure Pop1 and Pop2 are factors with levels in correct order
## Order pop1 and pop2 so they are the same order as that in allpops
## Set order pops should be plotted
allpops <- c("CLAC" ,"DUIN","LUIB","OBSE", "CLAM", "DUIM", "LUIM", "OBSM")

## Order pop1 and pop2 so they are the same order as that in allpops
adaptive_divergence_Fst <- apply(adaptive_divergence_Fst, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops))]
}
      ) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(adaptive_divergence_Fst)

## Convert Pop1 and Pop2 to factors with correct levels
adaptive_divergence_Fst <- adaptive_divergence_Fst %>%
  mutate(Pop1 = factor(V1, levels = allpops),
         Pop2 = factor(V2, levels = allpops))

## Triangle for outlier windows only
outlier_triangle <- data.frame(x = c(0.4, 8.495, 0.4, 0.4), y = c(0.505, 8.6, 8.6, 0.6))
## Plot heatmap of mean Fst between populations for outlier windows only
Fst_outlier_comp_plot <- ggplot() +
    geom_tile(data = adaptive_divergence_Fst[adaptive_divergence_Fst$anad_resi_outlier_99=="non-outlier",], aes(x = Pop1, y = Pop2, fill = mean_Fst)) +
    geom_tile(data = adaptive_divergence_Fst[adaptive_divergence_Fst$anad_resi_outlier_99=="outlier",], aes(x = Pop2, y = Pop1, fill = mean_Fst)) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90", name = "Fst", limits = c(0, 0.3)) +
    theme_classic() +
    geom_vline(xintercept = 4.5) +
    geom_hline(yintercept = 4.5) +
    geom_hline(yintercept = 8.5) +
    geom_vline(xintercept = 8.5) +
    geom_polygon(data = outlier_triangle, aes(x, y), fill = "NA", color = "grey30", linewidth = 1.5) +
    geom_polygon(data = outlier_triangle, aes(y, x), fill = "NA", color = "red", linewidth = 1.5) +
    theme(panel.background = element_rect(fill = "grey"), legend.frame = element_rect(colour = 'black'),
        axis.title = element_blank(), axis.text.y = element_text(angle=90, vjust = 4, hjust = 0.5),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        plot.margin = margin(b = 5.5, r = 5.5, t = 5.5, l = 30),
        axis.text.x = element_text(vjust=-1)) +
    annotate("text", x = -0.25, y = 6.5, label = "Migratory", angle = 90) + ## Axis annotation
    annotate("text", x = -0.25, y = 2.5, label = "Resident", angle = 90) +
    annotate("text", y = -0.25, x = 6.5, label = "Migratory") +
    annotate("text", y = -0.25, x = 2.5, label = "Resident") +
    annotate("text", y = 8.8, x = 2, label = "Non-outliers (99%)", col = "grey30") +
    annotate("text", y = 2, x = 8.8, label = "Outliers (1%)", angle = -90, col = "red") +
    coord_fixed(clip = 'off', x = c(0.5, 8.5), y = c(0.5, 8.5))

## Save plot
ggsave(hist_plot + Fst_outlier_comp_plot, file=paste0(fst_file, "_adaptive_divergence_Fst_comparison_99perc.png"), width=12, height=6)

