# Get vcf file from arguments
args <- commandArgs(trailingOnly=T)
geno.file <- args[1]
# geno.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/CLAC.genotypes.txt"
ROH.file <- args[2]
#  ROH.file <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/ROH/CLAC.out.RG.txt"

## load libraries
library(tidyverse)
library(data.table)
## Read in geno file
geno <- fread(geno.file, col.names = c("CHROM","POS","SAMPLE","GT"))
ROH  <- fread(ROH.file, skip=1,
              col.names = c("RG","Sample","Chromosome","Start","End",
                            "Length","Number_of_markers","Quality"))

# Get samples
individuals <- unique(geno$SAMPLE)

# Get combination of samples and ROH cutoffs
sample.cutoffs <- CJ(individuals, cutoff = seq(0, 1000000, 50000))

# Add interval columns for overlaps
geno[, `:=`(start = POS, end = POS)]
ROH[,  `:=`(start = Start, end = End)]

# Ensure keys for foverlaps
setkey(ROH, Chromosome, start, end)

# Iterate over combinations
ROH.calc <- rbindlist(lapply(1:nrow(sample.cutoffs), function(i) {
    samp   <- sample.cutoffs$individuals[i]
    cutoff <- sample.cutoffs$cutoff[i]

    # Subset ROH for this sample
    ROH.tmp <- ROH[Sample == samp & Length >= cutoff]
    nROH    <- nrow(ROH.tmp)
    sumROH  <- if (nROH > 0) sum(ROH.tmp$Length) else 0

    # Subset heterozygous genotypes
    geno.tmp <- geno[SAMPLE == samp & GT == "0/1"]
    
    nHet <- 0L
    if (nrow(ROH.tmp) > 0 && nrow(geno.tmp) > 0) {
        # Interval anti-join: SNPs NOT in ROH
        setkey(geno.tmp, CHROM, start, end)
        overlaps <- foverlaps(geno.tmp, ROH.tmp,
                              by.x = c("CHROM","start","end"),
                              by.y = c("Chromosome","start","end"),
                              nomatch = 0L, which = TRUE)
        geno.tmp <- geno.tmp[-overlaps$xid, ] # Remove overlaping regions
        nHet <- nrow(geno.tmp)
    } else {
        nHet <- nrow(geno.tmp)
    }

    data.table(samp, cutoff, nROH, sumROH, nHet)
}))

write.table(ROH.calc, file = paste0(gsub(".txt", "", ROH.file),".calcs.txt"), row.names = F, quote = F, col.names = TRUE, sep = ",")

