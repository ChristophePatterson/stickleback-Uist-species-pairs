# GO enrichment using KEGG following
# https://yulab-smu.top/biomedical-knowledge-mining-book/022-kegg.html
# https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html

## Install clusterProfiler
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationHub")
# BiocManager::install("ensembldb")
library(clusterProfiler)
library(AnnotationHub)
library(tidyverse)
library(patchwork)

## NEEDS TO BE RUN ON LAPTOP ##

# Colour Paletter
cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

## Three spined stickleback on KEGG
# https://www.kegg.jp/kegg-bin/show_organism?org=gat

# Setwd
setwd("/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Regions_of_interest")

sp_name <- "Gasterosteus aculeatus"
org_code <- "gat" # FOR KEGG

## Check that the organism is supported in KEGG
search_kegg_organism(org_code, by='kegg_code')

## Download annotation hub for specess
ah <- AnnotationHub()
ahs <- query(ah, sp_name)

ahs$genome

ahs
gas.acu <- ahs[['AH120714']] # Select list

# Perform enrichment analysis on test subset of gene
## Determine what the keytpe for that organism is usnig dummy gene
genelist <- c("120823408","120820219","120823426","120823570","120823573","120823472","120823563","120828214")
enrichKEGG(gene = genelist, organism = org_code, keyType = "kegg")
enrichGO(gene = genelist, OrgDb = gas.acu, keyType = "ENTREZID")

# Load in VEP output
vep.out <- read_table("stickleback_DUIN_minAC6_all_sig_top_regions.vep.out",comment = "#",
           col_names = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")) %>%
  mutate(Consequence= str_to_title(gsub("_", " ", Consequence)))


major.conseq.vars <- str_to_title(c("missense variant","inframe insertion",
  "inframe deletion","stop gained",
  "stop lost","frameshift variant",
  "start lost"))

vep.conseq <- vep.out %>%
  mutate(
    Consequence = factor(Consequence,levels = names(sort(table(vep.out$Consequence)))),
    Consequence.major = Consequence %in% major.conseq.vars) %>%
  dplyr::filter(Consequence.major) %>%
  group_by(Consequence) %>%
  summarise(n.SNPs = n(), n.genes = length(unique(Gene)))

p.conseq <- vep.out %>%
  mutate(
    Consequence = factor(Consequence,levels = names(sort(table(Consequence)))),
    Consequence.major = Consequence %in% major.conseq.vars) %>%
  dplyr::filter(Consequence.major) %>%
  ggplot() +
  geom_bar(aes(y = Consequence, fill = Consequence)) +
  geom_text(data = vep.conseq, aes(x = n.SNPs, y = Consequence, label = paste0(n.SNPs,"\n(", n.genes, ")")),
    hjust = -0.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)), name = "SNP Count") +
  scale_fill_manual(values = cbPalette) +
  theme_bw()


ggsave("Consequence_break_down.png", p.conseq)

# Read in annotation of DUKE genome
DUKE.annotated.gft <- read_table("GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly.gtf") %>%
  rename(c('name' = 'Gene'))

## How many genes are not identified
nrow(DUKE.annotated.gft)
length(na.omit(DUKE.annotated.gft$fGas.GeneID.match))/nrow(DUKE.annotated.gft)
length(na.omit(DUKE.annotated.gft$v5.GeneID.match))/nrow(DUKE.annotated.gft)
table(is.na(DUKE.annotated.gft$v5.GeneID.match)&is.na(DUKE.annotated.gft$fGas.GeneID.match))/nrow(DUKE.annotated.gft)

table(DUKE.annotated.gft$v5.GeneID.match==DUKE.annotated.gft$fGas.GeneID.match, useNA = "always")/nrow(DUKE.annotated.gft)


length(unique(DUKE.annotated.gft$fGas.name.match))

# Merge vep and DUke annotation
vep.out <- left_join(vep.out, DUKE.annotated.gft, by = "Gene")

# Read in all the genes that were succefully tetermined by blast from DUKE
DUKE.v5.successful.annotation <- read_table("Duke_GAcu_1.0_genomic_blast_v5GeneIDs.txt", col_names = "ID") %>%
  mutate(ID = as.character(ID))
DUKE.fGas.successful.annotation <- read_table("Duke_GAcu_1.0_genomic_blast_fGasGeneIDs.txt", col_names = "ID") %>%
  mutate(ID = as.character(ID))

# Get list of genes with major missense mutation filtering out NA's
fGas.GeneID.match.MajorConseq.uniq <- as.character(unique(na.omit(vep.out$fGas.GeneID.match[vep.out$Consequence %in% major.conseq.vars])))
v5.GeneID.match.MajorConseq.uniq <- as.character(unique(na.omit(vep.out$v5.GeneID.match[vep.out$Consequence %in% major.conseq.vars])))

# Run enrichment analysis for genes labeled from v5 and fGas
enrich.fGas <- enrichGO(gene = fGas.GeneID.match.MajorConseq.uniq , ont = "ALL",
                    keyType = "ENTREZID", OrgDb = gas.acu,  universe = DUKE.fGas.successful.annotation$ID)
enrich.v5 <- enrichGO(gene = fGas.GeneID.match.MajorConseq.uniq , ont = "ALL",
                    keyType = "ENTREZID", OrgDb = gas.acu,  universe = DUKE.v5.successful.annotation$ID) 

# Plot the GO terms
# For fGas
p <- barplot(enrich.fGas,  x = "Count",  split = "ONTOLOGY") +
  aes(fill = p.adjust) +
  facet_wrap(~ONTOLOGY, scales = "free",ncol =  1, space = "free_y",strip.position = "right") +
  scale_fill_gradient(name = "P-Adjusted", high = "skyblue", low = "firebrick", limits = c(0, 0.05)) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle("fGasAcu3.hap1.1")
  
p

# For v5
q <- barplot(enrich.v5, x = "Count",  split = "ONTOLOGY") +
  aes(fill = p.adjust) +
  facet_wrap(~ONTOLOGY, scales = "free",ncol =  1, space = "free_y", strip.position = "right") +
  scale_fill_gradient(name = "P-Adjusted", high = "skyblue", low = "firebrick", limits = c(0, 0.05)) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle("GAculeatus_UGA_version5")

q
# Combine
pq <- p + q + plot_layout(guides = "collect")

ggsave("GO_enrichment_analysis_fGas_and_v5.png", pq, width = 10, height = 12)

## Plot clusters of enriched genes
library(enrichplot)

cnet.fGas <- cnetplot(setReadable(enrich.fGas, gas.acu, 'ENTREZID'),
                      color_category = "skyblue",color_item = "orange",showCategory = 6)
cnet.v5 <- cnetplot(setReadable(enrich.v5, gas.acu, 'ENTREZID'), 
                    color_category = "skyblue",color_item = "orange",showCategory = 6) 

qcnet.fGas <- p + cnet.fGas
qcnet.v5 <- q + cnet.v5

ggsave("GO_enrichment_analysis_fGas_lv6cnet.png", qcnet.fGas, width = 15, height = 8)
ggsave("GO_enrichment_analysis_v5_lv6cnet.png", qcnet.v5, width = 15, height = 8)

cnet.fGas <- cnetplot(setReadable(enrich.fGas, gas.acu, 'ENTREZID'), 
                      color_category = "skyblue",color_item = "orange",
                      showCategory = 2)
cnet.v5 <- cnetplot(setReadable(enrich.v5, gas.acu, 'ENTREZID'), 
                    color_category = "skyblue",color_item = "orange",
                    showCategory = 2) 

qcnet.fGas <- p + cnet.fGas
qcnet.v5 <- q + cnet.v5

ggsave("GO_enrichment_analysis_fGas_lv2cnet.png", qcnet.fGas, width = 15, height = 8)
ggsave("GO_enrichment_analysis_v5_lv2cnet.png", qcnet.v5, width = 15, height = 8)


enriched.genes.fGas <-  str_split(enrich.fGas@result$geneID[enrich.fGas@result$p.adjust<=0.05],"/") %>%
  unlist() %>%
  unique()

enriched.genes.v5 <-  str_split(enrich.v5@result$geneID[enrich.v5@result$p.adjust<=0.05],"/") %>%
  unlist() %>%
  unique()

vep.out$enriched.fGas <- vep.out$fGas.GeneID.match %in% enriched.genes.fGas
vep.out$enriched.v5 <- vep.out$v5.GeneID.match %in% enriched.genes.v5

table(vep.out[vep.out$enriched.fGas&vep.out$Consequence=="missense_variant",]$fGas.GeneID.match)
table(vep.out[vep.out$enriched.fGas&vep.out$Consequence=="missense_variant",]$fGas.name.match)
table(vep.out[vep.out$enriched.v5&vep.out$Consequence=="missense_variant",]$v5.GeneID.match)
table(vep.out[vep.out$enriched.v5&vep.out$Consequence=="missense_variant",]$v5.name.match)
vep.out[vep.out$enriched.v5&vep.out$Consequence=="missense_variant",] %>%
  group_by(Gene) %>%
  summarise(chr = dplyr::first(chr), start = min(start), end = max(start), n.misssence = n()) %>%
  mutate(IGV.corr = paste0(chr, ":", start, "-", end))


# Run KEGG analysis
KEGG.fGas <- enrichKEGG(gene = fGas.GeneID.match.uniq, organism = org_code, keyType = "kegg")
KEGG.v5 <- enrichKEGG(gene = v5.GeneID.match.uniq, organism = org_code, keyType = "kegg")

head(KEGG.fGas)
head(KEGG.v5)
table(KEGG.fGas@result$GeneRatio)
KEGG.v5

vep.out.enriched <- vep.out[vep.out$fGas.GeneID.match%in%str_split(head(KEGG.fGas)$geneID, "/")[[1]],]
vep.out.enriched[vep.out.enriched$Consequence=="missense_variant",]

vep.out[grep("eda", vep.out$fGas.name.match),c("Location", "Consequence", "Gene", "fGas.name.match", "fGas.GeneID.match")]


KEGG.fGas@result[,c("ID", "geneID")]
browseKEGG(KEGG.v5,"gat00100")
browseKEGG(KEGG.v5,"gat04820")


