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
  mutate(gene_id = ifelse(type =="gene", name, 
                    ifelse(type =="transcript", str_split_i(name, "\\.", 1), gsub(";","",str_split_i(name, " ", 4)))),
         transcript_id = ifelse(type =="gene", NA, 
                            ifelse(type =="transcript", str_split_i(name, "\\.", 2), str_split_i(str_split_i(name, "\\.", 2), ";", 1))))

# Read in GasAcu3 annontation
fGas.gtf <- read_table("GCF_964276395.1/genomic.gff", comment = "#",col_names = FALSE) %>%
  rename(chr = X1, start = X4, end = X5, type = X3, name=X9) %>%
  filter(type == "mRNA") %>%
  mutate(gene.name = str_split_i(str_split_i(name, ";gene=", 2),";" ,1),
         parent = str_split_i(str_split_i(name, ";Parent=", 2),";" ,1),
         ID = str_split_i(str_split_i(name, "ID=", 2),";" ,1),
         GeneID = str_split_i(str_split_i(str_split_i(name, "GeneID:", 2),";" ,1),",",1))

# Read in v5 annontation
v5.gtf <- read_table("GCF_016920845.1/genomic.gff", comment = "#",col_names = FALSE) %>%
  rename(chr = X1, start = X4, end = X5, type = X3, name=X9) %>%
  filter(type == "mRNA") %>%
  mutate(#chr = factor(as.character(as.roman(as.numeric(gsub("chr", "", chr)))), levels = chr.order),
    gene.name = str_split_i(str_split_i(name, ";gene=", 2),";" ,1),
    parent = str_split_i(str_split_i(name, ";Parent=", 2),";" ,1),
    ID = str_split_i(str_split_i(name, "ID=", 2),";" ,1),
    GeneID = str_split_i(str_split_i(str_split_i(name, "GeneID:", 2),";" ,1),",",1)) 

## Matches with GasAcu3
DUKE_fGas_blast <- read_table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_vs_fGasAcu3_blastn_6.out", 
  col_names = c("transcript_id", "fGas_match", "per.match", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(gene = str_split_i(transcript_id, "\\.", 1),
         transcript = str_split_i(transcript_id, "\\.", 2))

# How many unique gene matches are for each gene and transcript?
DUKE_fGas_blast.unique <- DUKE_fGas_blast %>%
  group_by(gene, .add = TRUE) %>%
  group_by(transcript, .add = TRUE) %>%
  summarise(n = n_distinct(fGas_match),
            #max.per.match = max(per.match),
            #min.per.match = min(per.match),
            # mean.per.match = mean(per.match),
            best.bitscore = max(bitscore),
            best.evalue = min(evalue),
            best.hit.evalue = fGas_match[which.min(evalue)],
            best.hit.bitscore = fGas_match[which.max(bitscore)]) 
print(DUKE_fGas_blast.unique, width = Inf)

# How many unique gene matches are there for each gene?
DUKE_fGas_blast.unique.gene <- DUKE_fGas_blast.unique %>%
  group_by(gene) %>%
  summarise(n = n_distinct(best.hit.evalue),
            min.bitscore = min(best.bitscore),
            max.bitscore = max(best.bitscore),
            best.hit.bitscore = best.hit.bitscore[which.max(best.bitscore)]) %>%
  arrange(desc(n))

# Add the best blast match to the gtf file
DUKE.gtf$fGas_match_parent <- DUKE_fGas_blast.unique$best.hit.evalue[match(paste(DUKE.gtf$gene_id, DUKE.gtf$transcript_id, sep = "_"),
                                                                     paste(DUKE_fGas_blast.unique$gene, DUKE_fGas_blast.unique$transcript, sep = "_"))]

# For gene lines, match by gene_id only
DUKE.gtf$fGas_match_parent[DUKE.gtf$type == "gene"] <- DUKE_fGas_blast.unique.gene$best.hit.bitscore[match(DUKE.gtf$gene_id[DUKE.gtf$type == "gene"], DUKE_fGas_blast.unique.gene$gene)]
# How many genes do not have a match to fGasAcu3?
sum(is.na(DUKE.gtf$fGas_match_parent[DUKE.gtf$type == "gene"]))

## Matches with v5
DUKE_v5_blast <- read_table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_vs_v5_blastn_6.out", 
  col_names = c("transcript_id", "v5_match", "per.match", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(gene = str_split_i(transcript_id, "\\.", 1),
         transcript = str_split_i(transcript_id, "\\.", 2))

# DUKE_v5_blast$parent <- v5.gtf$parent[match(DUKE_v5_blast$v5_match, v5.gtf$ID)]

# How many unique gene matches are for each gene and transcript?
DUKE_v5_blast.unique <- DUKE_v5_blast %>%
  group_by(gene, .add = TRUE) %>%
  group_by(transcript, .add = TRUE) %>%
  summarise(n = n_distinct(v5_match),
            #max.per.match = max(per.match),
            #min.per.match = min(per.match),
            # mean.per.match = mean(per.match),
            best.bitscore = max(bitscore),
            best.evalue = min(evalue),
            best.hit.evalue = v5_match[which.min(evalue)],
            best.hit.bitscore = v5_match[which.max(bitscore)]) 
print(DUKE_v5_blast.unique, width = Inf)

# How many unique gene matches are there for each gene?
DUKE_v5_blast.unique.gene <- DUKE_v5_blast.unique %>%
  group_by(gene) %>%
  summarise(n = n_distinct(best.hit.evalue),
            min.bitscore = min(best.bitscore),
            max.bitscore = max(best.bitscore),
            best.hit.bitscore = best.hit.bitscore[which.max(best.bitscore)]) %>%
  arrange(desc(n))

# Add the best blast match to the gtf file
DUKE.gtf$v5_match_parent <- DUKE_v5_blast.unique$best.hit.evalue[match(paste(DUKE.gtf$gene_id, DUKE.gtf$transcript_id, sep = "_"),
                                                                     paste(DUKE_v5_blast.unique$gene, DUKE_v5_blast.unique$transcript, sep = "_"))]

# For gene lines, match by gene_id only
DUKE.gtf$v5_match_parent[DUKE.gtf$type == "gene"] <- DUKE_v5_blast.unique.gene$best.hit.bitscore[match(DUKE.gtf$gene_id[DUKE.gtf$type == "gene"], DUKE_v5_blast.unique.gene$gene)]
# How many genes do not have a match to v5?
sum(is.na(DUKE.gtf$v5_match_parent[DUKE.gtf$type == "gene"]))

## Match between the two blast outputs
sum(DUKE.gtf$fGas_match_parent == DUKE.gtf$v5_match_parent, na.rm = T)

# Add the gene name and gene ID from the fGas annotation to the gtf file based on the blast match
DUKE.gtf$fGas.name.match <- fGas.gtf$parent[match(DUKE.gtf$fGas_match_parent, gsub("rna-", "", fGas.gtf$ID))]
DUKE.gtf$fGas.GeneID.match <- fGas.gtf$GeneID[match(DUKE.gtf$fGas_match_parent, gsub("rna-", "", fGas.gtf$ID))]

# Add the gene name and gene ID from the v5 annotation to the gtf file based on the blast match
DUKE.gtf$v5.name.match <- v5.gtf$parent[match(DUKE.gtf$v5_match_parent, v5.gtf$ID)]
DUKE.gtf$v5.GeneID.match <- v5.gtf$GeneID[match(DUKE.gtf$v5_match_parent, v5.gtf$ID)]


# Write out the gtf file with the blast matches added as new columns
DUKE.gtf %>%
  select(chr, start, end, method, type, gene_id, transcript_id, fGas_match_parent, fGas.name.match, fGas.GeneID.match, v5_match_parent, v5.name.match, v5.GeneID.match) %>%
    mutate(fGas.name.match = gsub("gene-", "", fGas.name.match),
         v5.name.match = gsub("gene-", "", v5.name.match)) %>%
  write.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf", sep = "\t", quote = F, row.names = F, col.names = T)

# Subset to just the gene lines and write out a separate file
DUKE.gtf.genes <- DUKE.gtf %>%
  filter(type == "gene") %>%
  select(chr, start, end, name, fGas_match_parent, fGas.GeneID.match, fGas.name.match, v5_match_parent, v5.GeneID.match, v5.name.match) %>%
  mutate(fGas.name.match = gsub("gene-", "", fGas.name.match),
         v5.name.match = gsub("gene-", "", v5.name.match))

# Write out the gene-only gtf file
write.table(DUKE.gtf.genes, "GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly.gtf", sep = "\t", quote = F, row.names = F, col.names = T)

## How many gene matches are not the same between the two blast outputs?
sum(DUKE.gtf.genes$fGas.GeneID.match != DUKE.gtf.genes$v5.GeneID.match, na.rm = T) # 4841 
sum(DUKE.gtf.genes$fGas.GeneID.match != DUKE.gtf.genes$v5.GeneID.match, na.rm = T)/ nrow(DUKE.gtf.genes) # 0.22

# Print mismatches to a file
DUKE.gtf.genes %>%
  filter(fGas_match_parent != v5_match_parent) %>%
  write.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly_mismatches.gtf", sep = "\t", quote = F, row.names = F, col.names = T)

# How many genes are not ID to a name (have "LOC" in the name)
sum(grepl("LOC", DUKE.gtf.genes$fGas.name.match))
sum(grepl("LOC", DUKE.gtf.genes$v5.name.match))

# GeneIDs for Gene Ontology analysis
DUKE.gtf.genes %>%
  select(fGas.name.match) %>%
  filter(!is.na(fGas.name.match)) %>%
  write.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly_fGasGeneNames.txt", sep = "\t", quote = F, row.names = F, col.names = F)

DUKE.gtf.genes %>%
  select(v5.name.match) %>%
  filter(!is.na(v5.name.match)) %>%
  write.table("GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches_genesOnly_v5GeneNames.txt", sep = "\t", quote = F, row.names = F, col.names = F)