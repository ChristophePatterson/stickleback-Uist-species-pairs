# Rscriot tt.R <path to TT calculation> <output_folder>
## Code to produce custom plots from https://github.com/popgenDK/seqAfrica_giraffe
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly=T)

inTT <- paste0(args[1], "/TTcals/")
inSFS <- paste0(args[1], "/sfs/")
outprefix <- args[2]

pop_pair <- basename(outprefix)

TT.files <- list.files(inTT)
TTsfs.files <- list.files(inSFS)

TTcalc.files <- TT.files[grep(".years", TT.files)]

length(TTsfs.files)==length(TTcalc.files)

g.seg <- length(TTsfs.files)

TTcalcs <- do.call("rbind.data.frame", lapply(TTcalc.files, FUN = function(x) t(read.table(paste0(inTT,"/", x), nrows = 3))[2,]))
sfs.sum <- do.call("c", lapply(TTsfs.files, FUN = function(x) sum(read.table(paste0(inSFS,"/", x))[,2])))
TTcalcs <- cbind(TTcalc.files, TTcalcs, sfs.sum)
colnames(TTcalcs) <- c("TTparms","T1", "T2", "Na", "nsnps")


TTcalcs$range <- gsub(".tt.split.years", "", gsub(TTcalcs$TTparms, pattern = paste0(paste0(pop_pair, "_chr",1:22,"_"), collapse = '|'), replacement  = ""))
TTcalcs$chr <- stringr::str_split_i(TTcalcs$TTparms,"_", 3)
TTcalcs$chr <- factor(TTcalcs$chr, levels = unique(TTcalcs$chr)[order(as.numeric(gsub("chr", "", unique(TTcalcs$chr))))])
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
ggsave(paste0(outprefix, "_T1andT2_hist.png"), p+q, width = 10, height = 10)

p <- ggplot(TTcalcs) +
  geom_point(aes(x = start, y = T1/1000), col = "red") +
  geom_point(aes(x = start, y = T2/1000), col = "blue") +
  ylab("Year (ka)") +
  facet_grid(.~chr, scale = "free_x", space = "free_x") +
  ggtitle(pop_pair)

ggsave(paste0(outprefix, "_T1andT2.pdf"), p, width = 10, height = 5)
ggsave(paste0(outprefix, "_T1andT2.png"), p, width = 10, height = 5)


## Calculate weighted jackknife
weighted_jackknife <- function(values, weights) {
  ta <- values
  mj <- weights
  n_ <- sum(mj)
  ## Number of groups
  g <- length(ta)
  Oa <- mean(ta)
  # Calculate jacknight m(-j)
  a_mean <- 0
  i <- 1
  for(i in 1:g){
    a_count <- mj[i]
    a_term <- 1-1*(a_count/n_)
    a_mean <- a_mean+(a_term*ta[i])
  }
  # Calculate jackknife mean
  a_mean <- g*Oa-a_mean
  # Calculate jackknife var
  a_var <- 0
  for(i in 1:g){
    a_count <- mj[i]
    hj <- n_/a_count
    a_term <- hj*Oa-(hj-1)*ta[i]-a_mean
    a_var <- a_var+a_term/(hj-1)
  }
  a_var <- a_var/g
  
  return(c(obs_mean = Oa, jk_mean = a_mean, jk_var = a_var))
}


T1_JK <- weighted_jackknife(TTcalcs$T1, TTcalcs$nsnps)
T2_JK <- weighted_jackknife(TTcalcs$T2, TTcalcs$nsnps)

print(paste0(outprefix, "_T1_T2_jacknife_calc.txt"))

rbind(T1_JK,T2_JK)

write.table(rbind(T1_JK,T2_JK), paste0(outprefix, "_T1_T2_jacknife_calc.txt"))

