# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
geno.file <- args[1]
# geno.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/CLAC.genotypes.txt"
ROH.file <- args[2]
#  ROH.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/CLAC.out.RG.txt"

## load libraries
library(tidyverse)
## Readin geno.file
geno <- read_table(geno.file, col_names = c("CHROM","POS", "SAMPLE","GT"))
ROH <-  read_table(ROH.file, col_names = c("RG","Sample","Chromosome","Start","End","Length","Number of markers","Quality"), skip =1)

# Get samples
individuals <- unique(geno$SAMPLE)

#  Get combination of samples and ROH cutoffs
sample.cutoffs <- expand_grid(individuals, cutoff = seq(0, 1000000, 50000))

ROH.calc <- map2_dfr(sample.cutoffs$individuals, sample.cutoffs$cutoff, ~{
    samp <- .x
    cutoff <- .y
      
    ROH.tmp <- ROH[ROH$Sample==samp & ROH$Length>=cutoff,]
    # number of ROHs
    nROH <- nrow(ROH.tmp)
    # Length ROHs
    sumROH = sum(ROH.tmp$Length)
    # Het of nonROH regions
    geno.tmp <- geno %>%
                    filter(SAMPLE==samp & GT == "0/1" & CHROM %in% ROH.tmp$Chromosome) %>%
                    rowwise() %>%
                    filter(!any(CHROM == ROH.tmp$Chromosome &
                                POS >= ROH.tmp$Start &
                                POS <= ROH.tmp$End)) %>%
            ungroup()

    nHet = nrow(geno.tmp)
    return(data.frame(samp, cutoff, nROH, sumROH, nHet))
})

write.table(ROH.calc, file = paste0(gsub(".txt", "" ROH.file),".calcs.txt", row.names = F, quote = F, col.names = TRUE, sep = ","))