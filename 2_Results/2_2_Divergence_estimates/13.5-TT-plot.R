# Rscriot tt.R <path to TT calculation> <output_folder>
## Taken from https://github.com/popgenDK/seqAfrica_giraffe
library(tidyverse)

args <- commandArgs(trailingOnly=T)

inTT <- args[1]
outprefix <- args[2]

TTcalcs.files <- list.files(inTT)

TTcalcs.files <- TTcalcs.files[grep(".years", TTcalcs.files)]

TTcalcs <- do.call("rbind.data.frame", lapply(TTcalcs.files, FUN = function(x) t(read.table(paste0(inTT,"/", x), nrows = 3))[2,]))

TTcalcs <- cbind(TTcalcs.files, TTcalcs)

colnames(TTcalcs) <- c("TTparms","T1", "T2", "Na")

TTcalcs$range <- gsub(".tt.split.years", "",gsub("OBSE_OBSM_chr1_", "", TTcalcs$TTparms))

TTcalcs$start <- as.numeric(stringr::str_split_i(TTcalcs$range,"-", 1))
TTcalcs$end <- as.numeric(stringr::str_split_i(TTcalcs$range,"-", 2))
TTcalcs$T1 <- as.numeric(TTcalcs$T1)
TTcalcs$T2 <- as.numeric(TTcalcs$T2)

## Reorder calcules
TTcalcs <- TTcalcs[order(TTcalcs$start),]

p <- ggplot(TTcalcs) +
    geom_histogram(aes(x = T1), fill = "blue", alpha = 0.5) +
    geom_histogram(aes(x = T2), fill = "red", alpha = 0.5) +
    ggtitle(outprefix)

ggsave(paste0(outprefix, "_T1andT2_hist.pdf"), p, width = 10, height = 10)

p <- ggplot(TTcalcs) +
    geom_point(aes(x = start, y = T1/1000), col = "red") +
    geom_point(aes(x = start, y = T2/1000), col = "blue") +
    ylab("Year (ka)") +
    ggtitle(outprefix)

ggsave(paste0(outprefix, "_T1andT2.pdf"), p, width = 10, height = 5)

