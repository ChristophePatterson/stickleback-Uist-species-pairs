library("tidyverse")

setwd("/gpfs01/home/mbzcp2/data/sticklebacks/genomes")

## Read in chrom info
chr <- as_tibble(read.table("GCA_046562415.1/GCA_046562415.1_seq_report.tsv", sep = "\t", header = T))
chr$Sequence.name <- gsub("chr", "", chr$Sequence.name)
# Calculate cummulative sum
chr$Cum.Seq.length <- c(0, cumsum(chr$Seq.length[1:(nrow(chr)-1)]))
chr.order <- chr$Sequence.name

# Read in DUKE annotation
DUKE.gtf <- read.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_ChrNames_fixed_notab.gtf", sep = "\t") %>%
  tibble() %>%
  rename(chr = V1, start = V4, end = V5, method = V2, type = V3, name = V9) %>%
  # filter(type == "gene") %>%
  mutate(transcript_id = gsub(";","",str_split_i(name, " ", 2)),
         gene_id = gsub(";","",str_split_i(name, " ", 4)))

# Read in GasAcu3 annontation
fGas.gtf <- read_table("GCF_964276395.1/genomic.gff", comment = "#",col_names = FALSE) %>%
  rename(chr = X1, start = X4, end = X5, type = X3, name=X9) %>%
  mutate(gene.name = str_split_i(str_split_i(name, ";gene=", 2),";" ,1),
         parent = str_split_i(str_split_i(name, ";Parent=", 2),";" ,1),
         ID = str_split_i(str_split_i(name, "ID=", 2),";" ,1)) 

# Read in v5 annontation
v5.gtf <- read_table("GCF_016920845.1/genomic.gff", comment = "#",col_names = FALSE) %>%
  rename(chr = X1, start = X4, end = X5, type = X3, name=X9) %>%
  mutate(#chr = factor(as.character(as.roman(as.numeric(gsub("chr", "", chr)))), levels = chr.order),
    gene.name = str_split_i(str_split_i(name, ";gene=", 2),";" ,1),
    parent = str_split_i(str_split_i(name, ";Parent=", 2),";" ,1),
    ID = str_split_i(str_split_i(name, "ID=", 2),";" ,1)) 

## Matches with GasAcu3
DUKE_fGas_blast <- read.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_vs_fGasAcu3_blastn.out", sep = "\t") %>%
  tibble() %>%
  rename(transcript_id = V1, fGas_match = V2, per.match = V3) %>%
  mutate(gene = str_split_i(transcript_id, "\\.", 1))

DUKE_fGas_blast$parent <- fGas.gtf$parent[match(DUKE_fGas_blast$fGas_match, gsub("rna-", "",fGas.gtf$ID))]

DUKE.gtf$fGas_match_parent <- NA
DUKE.gtf$fGas_match_parent[DUKE.gtf$type=="transcript"] <- (DUKE_fGas_blast$parent[match(DUKE.gtf$name, DUKE_fGas_blast$transcript_id)])[DUKE.gtf$type=="transcript"]
DUKE.gtf$fGas_match_parent[DUKE.gtf$type=="gene"] <- (DUKE_fGas_blast$parent[match(DUKE.gtf$name, DUKE_fGas_blast$gene)])[DUKE.gtf$type=="gene"]
DUKE.gtf$fGas_match_parent[!DUKE.gtf$type%in%c("gene","transcript")] <- (DUKE_fGas_blast$parent[match(DUKE.gtf$transcript_id, DUKE_fGas_blast$transcript_id)])[!DUKE.gtf$type%in%c("gene","transcript")]

## Matches with v5
DUKE_v5_blast <- read.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_vs_v5_blastn.out", sep = "\t") %>%
  tibble() %>%
  rename(transcript_id = V1, v5_match = V2, per.match = V3) %>%
  mutate(gene = str_split_i(transcript_id, "\\.", 1))

DUKE_v5_blast$parent <- v5.gtf$parent[match(DUKE_v5_blast$v5_match, v5.gtf$ID)]

DUKE.gtf$v5_match_parent <- NA
DUKE.gtf$v5_match_parent[DUKE.gtf$type=="transcript"] <- (DUKE_v5_blast$parent[match(DUKE.gtf$name, DUKE_v5_blast$transcript_id)])[DUKE.gtf$type=="transcript"]
DUKE.gtf$v5_match_parent[DUKE.gtf$type=="gene"] <- (DUKE_v5_blast$parent[match(DUKE.gtf$name, DUKE_v5_blast$gene)])[DUKE.gtf$type=="gene"]
DUKE.gtf$v5_match_parent[!DUKE.gtf$type%in%c("gene","transcript")] <- (DUKE_v5_blast$parent[match(DUKE.gtf$transcript_id, DUKE_v5_blast$transcript_id)])[!DUKE.gtf$type%in%c("gene","transcript")]
 

# Write out the gtf file with the blast matches added as new columns
write.table(DUKE.gtf, "GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf", sep = "\t", quote = F, row.names = F, col.names = T)

# Subset to just the gene lines and write out a separate file
DUKE.gtf.genes <- DUKE.gtf %>%
  filter(type == "gene") %>%
  select(chr, start, end, name, fGas_match_parent, v5_match_parent) %>%
  mutate(fGas_match_parent = gsub("gene-", "", fGas_match_parent),
         v5_match_parent = gsub("gene-", "", v5_match_parent))

write.table(DUKE.gtf.genes, "GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly.gtf", sep = "\t", quote = F, row.names = F, col.names = T)

## How many gene matches are not the same between the two blast outputs?
sum(DUKE.gtf.genes$fGas_match_parent != DUKE.gtf.genes$v5_match_parent, na.rm = T) # 4841 
