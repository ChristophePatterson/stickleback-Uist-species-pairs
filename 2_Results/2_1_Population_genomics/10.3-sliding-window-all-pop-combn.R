library(tidyverse)

sliding_input <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/sliding-window/"

### ## CLAC resi vs all non CLAC anad
### CLACr_v_nCLACa <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_allanad.csv"))
### colnames(CLACr_v_nCLACa)[grep("Fst",colnames(CLACr_v_nCLACa))] <- "Fst"
### 
### ## nonCLAC resi vs nonCLAC anad
### nCLACr_v_nCLAMa <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_nCLAC_resi_nCLAManad.csv"))
### colnames(nCLACr_v_nCLAMa)[grep("Fst",colnames(nCLACr_v_nCLAMa))] <- "Fst"
### 
### CLACr_v_nCLACr <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_nCLACresi.csv"))
### colnames(CLACr_v_nCLACr)[grep("Fst",colnames(CLACr_v_nCLACr))] <- "Fst"
### 
scaf <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
### CLACr_v_nCLACa$chr <- scaf$Sequence.name[match(CLACr_v_nCLACa$scaffold,scaf$RefSeq.seq.accession)]
### nCLACr_v_nCLAMa$chr <- scaf$Sequence.name[match(nCLACr_v_nCLAMa$scaffold,scaf$RefSeq.seq.accession)]
### CLACr_v_nCLACr$chr <- scaf$Sequence.name[match(CLACr_v_nCLACr$scaffold,scaf$RefSeq.seq.accession)]
### 
### 
### p <- ggplot() +
###   geom_line(data = CLACr_v_nCLACa, aes(x = mid, y = Fst, group = chr, col = "CLACr_v_nCLAMa")) +
###   geom_line(data = nCLACr_v_nCLAMa, aes(x = mid, y = -Fst, group = chr, col = "nCLACr_v_nCLAMa")) +
###   scale_color_manual(values = c("skyblue", "orange")) +
###   facet_grid(chr~., scale = "free_x", space = "free_x") +
###   theme_classic() +
###   theme(panel.spacing = unit(0,'lines'), legend.position = "bottom")
### 
### ggsave(paste0("test.pdf"), p, width = 10, height = 30)
### ggsave(paste0(sliding_input,"Mirrored_CLACr_vs_nonCLACa_Fst.pdf"), p, width = 10, height = 30)
### 
### p <- ggplot() +
###   geom_line(data = CLACr_v_nCLACa, aes(x = mid, y = Fst, group = chr, col = "CLAC_resi_v_non-CLAM_anad"),alpha = 0.8) +
###   geom_line(data = nCLACr_v_nCLAMa, aes(x = mid, y = -Fst, group = chr, col = "non-CLAC_resi_v_non-CLAM-anad"),alpha = 0.8) +
###   geom_line(data = CLACr_v_nCLACr, aes(x = mid, y = -Fst, group = chr, col = "CLAC_resi_v_non-CLAC_resi"),alpha = 0.8) +
###   scale_color_manual(values = c("skyblue", "orange", "darkred")) +
###   facet_grid(chr~., scale = "free_x", space = "free_x") +
###   theme_classic() +
###   theme(panel.spacing = unit(0,'lines'), legend.position = "top")
### 
### ggsave(paste0("test.pdf"), p, width = 10, height = 30)
### ggsave(paste0(sliding_input,"Mirrored_CLACr_vs_nonCLAMa_vs_nCLACr.pdf"), p, width = 10, height = 30)

## CLAC_CLAM
CLAC_CLAM <- read.csv(paste0(sliding_input,"/All_Pop_comparison/CLAC_CLAM/sliding_window_w25kb_s5kb_m1_CLAC_CLAM.csv"))
colnames(CLAC_CLAM)[grep("Fst",colnames(CLAC_CLAM))] <- "Fst"
colnames(CLAC_CLAM)[grep("dxy",colnames(CLAC_CLAM))] <- "dxy"
colnames(CLAC_CLAM)[grep("pi_CLAC",colnames(CLAC_CLAM))] <- "pi_resi"
colnames(CLAC_CLAM)[grep("pi_CLAM",colnames(CLAC_CLAM))] <- "pi_anad"
CLAC_CLAM$chr <- scaf$Sequence.name[match(CLAC_CLAM$scaffold,scaf$RefSeq.seq.accession)]
CLAC_CLAM$pops <- "CLAC_CLAM"

## DUIN_DUIM
DUIN_DUIM <- read.csv(paste0(sliding_input,"/All_Pop_comparison/DUIM_DUIN/sliding_window_w25kb_s5kb_m1_DUIM_DUIN.csv"))
colnames(DUIN_DUIM)[grep("Fst",colnames(DUIN_DUIM))] <- "Fst"
colnames(DUIN_DUIM)[grep("dxy",colnames(DUIN_DUIM))] <- "dxy"
colnames(DUIN_DUIM)[grep("pi_DUIN",colnames(DUIN_DUIM))] <- "pi_resi"
colnames(DUIN_DUIM)[grep("pi_DUIM",colnames(DUIN_DUIM))] <- "pi_anad"
DUIN_DUIM$chr <- scaf$Sequence.name[match(DUIN_DUIM$scaffold,scaf$RefSeq.seq.accession)]
DUIN_DUIM$pops <- "DUIN_DUIM"

## OBSE_OBSM
OBSE_OBSM <- read.csv(paste0(sliding_input,"/All_Pop_comparison/OBSE_OBSM/sliding_window_w25kb_s5kb_m1_OBSE_OBSM.csv"))
colnames(OBSE_OBSM)[grep("Fst",colnames(OBSE_OBSM))] <- "Fst"
colnames(OBSE_OBSM)[grep("dxy",colnames(OBSE_OBSM))] <- "dxy"
colnames(OBSE_OBSM)[grep("pi_OBSE",colnames(OBSE_OBSM))] <- "pi_resi"
colnames(OBSE_OBSM)[grep("pi_OBSM",colnames(OBSE_OBSM))] <- "pi_anad"
OBSE_OBSM$chr <- scaf$Sequence.name[match(OBSE_OBSM$scaffold,scaf$RefSeq.seq.accession)]
OBSE_OBSM$pops <- "OBSE_OBSM"

## LUIB_LUIM
LUIB_LUIM <- read.csv(paste0(sliding_input,"/All_Pop_comparison/LUIB_LUIM/sliding_window_w25kb_s5kb_m1_LUIB_LUIM.csv"))
colnames(LUIB_LUIM)[grep("Fst",colnames(LUIB_LUIM))] <- "Fst"
colnames(LUIB_LUIM)[grep("dxy",colnames(LUIB_LUIM))] <- "dxy"
colnames(LUIB_LUIM)[grep("pi_LUIB",colnames(LUIB_LUIM))] <- "pi_resi"
colnames(LUIB_LUIM)[grep("pi_LUIM",colnames(LUIB_LUIM))] <- "pi_anad"
LUIB_LUIM$chr <- scaf$Sequence.name[match(LUIB_LUIM$scaffold,scaf$RefSeq.seq.accession)]
LUIB_LUIM$pops <- "LUIB_LUIM"

Pop_Fst <- rbind(OBSE_OBSM, CLAC_CLAM, LUIB_LUIM, DUIN_DUIM)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = Fst, group = chr, col = pops)) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$Fst))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave(paste0(sliding_input,"/Fst_all_pop_pairs.pdf"), p, width = 10, height = 30)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = dxy, group = chr, col = pops), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$dxy))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave(paste0(sliding_input,"/dxy_all_pop_pairs.pdf"), p, width = 10, height = 30)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = pi_resi, group = chr, col = pops), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$pi_resi))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave(paste0(sliding_input,"/pi_resi_all_pop_pair_comparisons.pdf"), p, width = 10, height = 30)


p <- ggplot(Pop_Fst[sample(nrow(Pop_Fst)),]) +
  geom_point(aes(x = dxy, y = Fst, col = pops), alpha = 0.5) +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave(paste0(sliding_input,"/Fst_Dxy_pop_pair_comparisons.pdf"), p, width = 10, height = 10)
ggsave(paste0(sliding_input,"/Fst_Dxy_pop_pair_comparisons.png"), p, width = 10, height = 10)

