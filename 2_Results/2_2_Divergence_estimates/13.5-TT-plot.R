# Rscriot tt.R <path to TT calculation> <output_folder>
## Code to produce custom plots from https://github.com/popgenDK/seqAfrica_giraffe
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly=T)

inTT <- args[1]
outprefix <- args[2]

pop_pair <- basename(outprefix)

TTcalcs.files <- list.files(inTT)

TTcalcs.files <- TTcalcs.files[grep(".years", TTcalcs.files)]

TTcalcs <- do.call("rbind.data.frame", lapply(TTcalcs.files, FUN = function(x) t(read.table(paste0(inTT,"/", x), nrows = 3))[2,]))

TTcalcs <- cbind(TTcalcs.files, TTcalcs)

colnames(TTcalcs) <- c("TTparms","T1", "T2", "Na")

head(TTcalcs$TTparms)

TTcalcs$range <- gsub(".tt.split.years", "", gsub(TTcalcs$TTparms, pattern = paste0(paste0(pop_pair, "_chr",1:22,"_"), collapse = '|'), replacement  = ""))
TTcalcs$chr <- stringr::str_split_i(TTcalcs$TTparms,"_", 3)

TTcalcs$start <- as.numeric(stringr::str_split_i(TTcalcs$range,"-", 1))
TTcalcs$end <- as.numeric(stringr::str_split_i(TTcalcs$range,"-", 2))
TTcalcs$T1 <- as.numeric(TTcalcs$T1)
TTcalcs$T2 <- as.numeric(TTcalcs$T2)

## Reorder calcules
TTcalcs <- TTcalcs[order(TTcalcs$start),]

summary(TTcalcs$T1)
summary(TTcalcs$T2)

p <- ggplot(TTcalcs) +
    geom_histogram(aes(x = T1), fill = "blue", alpha = 0.5) +
    geom_histogram(aes(x = T2), fill = "red", alpha = 0.5) +
    ggtitle(pop_pair)

q <- ggplot(TTcalcs) +
    geom_boxplot(aes(y = T1, x = "T1"), outlier.shape = NA) +
    geom_boxplot(aes(y = T2, x = "T2"), outlier.shape = NA) +
    ggtitle(pop_pair)   

ggsave(paste0(outprefix, "_T1andT2_hist.pdf"), p+q, width = 10, height = 10)

p <- ggplot(TTcalcs) +
    geom_point(aes(x = start, y = T1/1000), col = "red") +
    geom_point(aes(x = start, y = T2/1000), col = "blue") +
    ylab("Year (ka)") +
    facet_grid(.~chr, scale = "free_x", space = "free_x") +
    ggtitle(pop_pair)

ggsave(paste0(outprefix, "_T1andT2.pdf"), p, width = 10, height = 5)

