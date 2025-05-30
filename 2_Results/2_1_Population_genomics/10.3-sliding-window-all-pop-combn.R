library(tidyverse)
library(vcfR)

sliding_input <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/sliding-window/"

## CLAC resi vs all non CLAC anad
CLACr_v_nCLACa <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_allanad.csv"))
colnames(CLACr_v_nCLACa)[grep("Fst",colnames(CLACr_v_nCLACa))] <- "Fst"

## nonCLAC resi vs nonCLAC anad
nCLACr_v_nCLAMa <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_nCLAC_resi_nCLAManad.csv"))
colnames(nCLACr_v_nCLAMa)[grep("Fst",colnames(nCLACr_v_nCLAMa))] <- "Fst"

CLACr_v_nCLACr <- read.csv(paste0(sliding_input,"/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_nCLACresi.csv"))
colnames(CLACr_v_nCLACr)[grep("Fst",colnames(CLACr_v_nCLACr))] <- "Fst"

scaf <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
CLACr_v_nCLACa$chr <- scaf$Sequence.name[match(CLACr_v_nCLACa$scaffold,scaf$RefSeq.seq.accession)]
nCLACr_v_nCLAMa$chr <- scaf$Sequence.name[match(nCLACr_v_nCLAMa$scaffold,scaf$RefSeq.seq.accession)]
CLACr_v_nCLACr$chr <- scaf$Sequence.name[match(CLACr_v_nCLACr$scaffold,scaf$RefSeq.seq.accession)]


p <- ggplot() +
  geom_line(data = CLACr_v_nCLACa, aes(x = mid, y = Fst, group = chr, col = "CLACr_v_nCLAMa")) +
  geom_line(data = nCLACr_v_nCLAMa, aes(x = mid, y = -Fst, group = chr, col = "nCLACr_v_nCLAMa")) +
  scale_color_manual(values = c("skyblue", "orange")) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "bottom")

ggsave(paste0("test.pdf"), p, width = 10, height = 30)
ggsave(paste0(sliding_input,"Mirrored_CLACr_vs_nonCLACa_Fst.pdf"), p, width = 10, height = 30)

p <- ggplot() +
  geom_line(data = CLACr_v_nCLACa, aes(x = mid, y = Fst, group = chr, col = "CLAC_resi_v_non-CLAM_anad"),alpha = 0.8) +
  geom_line(data = nCLACr_v_nCLAMa, aes(x = mid, y = -Fst, group = chr, col = "non-CLAC_resi_v_non-CLAM-anad"),alpha = 0.8) +
  geom_line(data = CLACr_v_nCLACr, aes(x = mid, y = -Fst, group = chr, col = "CLAC_resi_v_non-CLAC_resi"),alpha = 0.8) +
  scale_color_manual(values = c("skyblue", "orange", "darkred")) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave(paste0("test.pdf"), p, width = 10, height = 30)
ggsave(paste0(sliding_input,"Mirrored_CLACr_vs_nonCLAMa_vs_nCLACr.pdf"), p, width = 10, height = 30)

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

Pop_Fst <- rbind(CLAC_CLAM, DUIN_DUIM, OBSE_OBSM, LUIB_LUIM)

p <- ggplot() +
  geom_line(data = CLAC_CLAM, aes(x = mid, y = Fst, group = chr, col = "CLAC_CLAM")) +
  geom_line(data = LUIB_LUIM, aes(x = mid, y = -Fst, group = chr, col = "LUIB_LUIM")) +
  scale_color_manual(values = c("skyblue", "orange")) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave("test.pdf", p, width = 10, height = 30)

Pop_Fst <- rbind(OBSE_OBSM, CLAC_CLAM, LUIB_LUIM)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = Fst, group = chr, col = pops), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$Fst))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave("test.pdf", p, width = 10, height = 30)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = dxy, group = chr, col = pops), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$dxy))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave("dxy_all_pop_comparisons.pdf", p, width = 10, height = 30)

p <- ggplot(Pop_Fst) +
  geom_line(aes(x = mid, y = pi_resi, group = chr, col = pops), alpha = 0.5) +
  facet_grid(chr~., scale = "free_x", space = "free_x") +
  theme_classic() +
  ylim(c(-0.1, max(Pop_Fst$pi_resi))) +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave("pi_resi_all_pop_comparisons.pdf", p, width = 10, height = 30)

##### Calculation of private alleles #####

vcf.SNPs <- read.vcfR("/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz")

## Calculated sequecning error rate from duplicated samples
technical_dups <- vcf.SNPs[samples = c("Obsm_640", "Obsm_641")]

technical_dups_gt <- extract.gt(technical_dups)
sum(as.numeric(na.omit(!technical_dups_gt[,1]==technical_dups_gt[,2])))/length(na.omit(technical_dups_gt[,1]==technical_dups_gt[,2]))*100

technical_dups_gt_errors <- na.omit(technical_dups_gt[technical_dups_gt[,1]!=technical_dups_gt[,2],])
table(paste0(technical_dups_gt_errors[,1], "-", technical_dups_gt_errors[,2]))

## Remove one of the dupicalted samples
vcf.SNPs <- vcf.SNPs[samples = colnames(vcf.SNPs@gt)[colnames(vcf.SNPs@gt)!="Obsm_641"]]
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs)]

## Read in sequence data populations infor
seq_data <- read.csv("/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv", header = T)
seq_data <- seq_data[match(colnames(vcf.SNPs@gt)[-1], seq_data$individual),]
seq_data$individual==(colnames(vcf.SNPs@gt)[-1])
seq_data$Population[seq_data$Region=="Maine"] <- "Maine"

# List of populations that I am interested in
Uist_pops <- c("OBSE", "OBSM", "OLAM", "OLAV", "OLST", "CLAM", "CLAC", "LUIB", "LUIM", "DUIM", "DUIN", "TORM", "TOST")
## Subset just to UIST data
seq_data_Uist <- seq_data[seq_data$Population%in%Uist_pops,]

## Subset vcf just to Uist
vcf.SNPs.Uist <- vcf.SNPs[samples = seq_data_Uist$individual]
vcf.SNPs.Uist <- vcf.SNPs.Uist[is.biallelic(vcf.SNPs.Uist)]
vcf.SNPs.Uist <- vcf.SNPs.Uist[is.polymorphic(vcf.SNPs.Uist, na.omit = T)]

# Remove non-CLAC samples
pop <- "CLAC"
vcf.tmp <- vcf.SNPs.Uist[samples = seq_data_Uist$individual[seq_data_Uist$Population!="CLAC"]]
# Retain samples that are no longer variable (so private to CLAC)
vcf.tmp <- vcf.tmp[!is.polymorphic(vcf.tmp,na.omit = T)]
## Poistion of these private allelle
priv_allele_pos <- paste(vcf.tmp@fix[,1], vcf.tmp@fix[,2])
## Number of private alleles
dim(vcf.tmp[1])
## Create vcf of CLACs private alleles
# Subset just to clac
vcf.CLAC.priv <- vcf.SNPs.Uist[samples = seq_data_Uist$individual[seq_data_Uist$Population=="CLAC"]]
## All snps
all.pos <- paste(vcf.CLAC.priv@fix[,1], vcf.CLAC.priv@fix[,2])
## Remove snps that arn't private
vcf.CLAC.priv <- vcf.CLAC.priv[(all.pos%in%priv_allele_pos),]
any(!vcf.tmp@fix[,2]==vcf.CLAC.priv@fix[,2])

pos.CLAC.priv <- as.data.frame(vcf.CLAC.priv@fix[,c(1,2)])
scaf <- read.table("/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv", sep = "\t", header = T)
pos.CLAC.priv$chr <- scaf$Sequence.name[match(pos.CLAC.priv$CHROM,scaf$RefSeq.seq.accession)]

p <- ggplot(pos.CLAC.priv) +
  geom_jitter(aes(as.numeric(POS),chr), width = 0, height = 0.1)

p <- ggplot(pos.CLAC.priv,aes(as.numeric(POS))) +
  geom_line(data = CLAC_CLAM, aes(x = mid, y = Fst*60, group = chr, col = "Fst CLAC_CLAM"), alpha = 0.5) +
  geom_line(data = CLACr_v_nCLACa, aes(x = mid, y = -Fst*60, group = chr, col = "Fst CLACr_v_nCLAMa"), alpha = 0.5) +
  geom_line(data = CLACr_v_nCLACr, aes(x = mid, y = -Fst*60, group = chr, col = "Fst CLACr_v_nCLAMr"), alpha = 0.5) +
  geom_freqpoly(aes(as.numeric(POS), col = "Freq Private Alleles CLAC"), binwidth = 25000, alpha = 0.5) +
  scale_y_continuous(name = "First Axis",sec.axis = sec_axis( trans=~./60, name="Second Axis")) +
  scale_color_manual(values = c("skyblue", "deepskyblue", "orange", "darkred")) +
  facet_grid(chr~.) +
  theme_classic() +
  theme(panel.spacing = unit(0,'lines'), legend.position = "top")

ggsave("test.pdf", p , width = 10, height =  30)


## Loop that calculate number of private alleles for each population
priv_alleles <- list()
for(pop in unique(seq_data_Uist$Population)){
  vcf.tmp <- vcf.SNPs.Uist[samples = seq_data_Uist$individual[seq_data_Uist$Population!=pop]]
  vcf.tmp <- vcf.tmp[!is.polymorphic(vcf.tmp,na.omit = T)]
  priv_alleles[[pop]] <- c(pops = pop, dim(vcf.tmp))
  print(pop)
}

priv_alleles_df <- as.data.frame(do.call("rbind",priv_alleles))
