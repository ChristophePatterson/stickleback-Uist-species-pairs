# Rscript 10.5-Private-alleles-plot.R <path to bcfstats>
## Code to produce custom plots from https://github.com/popgenDK/seqAfrica_giraffe
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly=T)

popPath <- args[1]

inAllStats <- paste0(popPath, ".PSC.stats.txt")
inPrivStats <- paste0(popPath, ".PSC.priv.stats.txt")

inAllStats <- read.table(inAllStats, skip = 2)
colnames(inAllStats) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions","nIndels",
                            "average_depth","nSingletons","nHapRef","nHapAlt","nMissing")

inPrivStats <- read.table(inPrivStats, skip = 2)
colnames(inPrivStats) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions","nIndels",
                            "average_depth","nSingletons","nHapRef","nHapAlt","nMissing")


inAllStats$pHets <- inAllStats$nHets/(inAllStats$nRefHom+inAllStats$nNonRefHom)
inPrivStats$pHets <- inPrivStats$nHets/(inPrivStats$nRefHom+inPrivStats$nNonRefHom)

p <- ggplot() +
        geom_boxplot(data = inAllStats, aes("All SNPs", pHets), outlier.shape =  NA) +
        geom_jitter(data = inAllStats, aes("All SNPs", pHets), height =  0) +
        geom_boxplot(data = inPrivStats, aes("Private SNPs", pHets), outlier.shape =  NA) +
        geom_jitter(data = inPrivStats, aes("Private SNPs", pHets), height =  0) +
        ggtitle(basename(popPath))

ggsave(paste0(popPath,"_het_all_v_priv.pdf"), p)
