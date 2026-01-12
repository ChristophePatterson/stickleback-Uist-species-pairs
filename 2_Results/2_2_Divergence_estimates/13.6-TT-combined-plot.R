library(ggplot2)
args <- commandArgs(trailingOnly=T)

inTT <- paste0(args[1], "/plots")

TT.files <- list.files(inTT)

TTcalc.files <- TT.files[grep(".txt", TT.files)]



TTcalcs <- do.call("rbind", lapply(TTcalc.files, FUN = function(x) read.table(paste0(inTT,"/", x))))
TTcalcs$T1_2 <- c("T1", "T2")
TTcalcs$Pops <- gsub("_T1_T2_jacknife_calc.txt", "", basename(TTcalc.files))
TTcalcs$Pop1 <- stringr::str_split_i(TTcalcs$Pops,"_", 1)
TTcalcs$Pop2 <- stringr::str_split_i(TTcalcs$Pops,"_", 2)
TTcalcs$Ecotype1 <- NA
TTcalcs$Ecotype2 <- NA
TTcalcs$Ecotype1[TTcalcs$Pop1%in%c("LUIM","CLAM","DUIM","OBSM","OLAM")] <- "anad"
TTcalcs$Ecotype2[TTcalcs$Pop2%in%c("LUIM","CLAM","DUIM","OBSM","OLAM")] <- "anad"
TTcalcs$Ecotype1[TTcalcs$Pop1%in%c("LUIB","CLAC","DUIN","OBSE")] <- "resi"
TTcalcs$Ecotype2[TTcalcs$Pop2%in%c("LUIB","CLAC","DUIN","OBSE")] <- "resi"
TTcalcs$Ecotype1[TTcalcs$Pop1%in%c("OLAV")] <- "fw"
TTcalcs$Ecotype2[TTcalcs$Pop2%in%c("OLAV")] <- "fw"


TTcalcs
TTcalcs$Ecotype_pair <- apply(TTcalcs[,c("Ecotype1", "Ecotype2")], MARGIN = 1, FUN = function(x) paste(sort(c(x[1], x[2])), collapse = "-"))

unique(TTcalcs$Ecotype_pair)

p <- ggplot(TTcalcs) +
    geom_boxplot(aes(Ecotype_pair, jk_mean, col = Ecotype_pair), outlier.shape = NA) +
    geom_jitter(aes(Ecotype_pair, jk_mean, col = Ecotype_pair), height = 0, width = 0.25) +
    facet_wrap(~T1_2)

ggsave(file = "test.pdf", p, width = 10, height = 7)
ggsave(file = "test.png", p, width = 10, height = 7)

q <- ggplot(TTcalcs) +
    geom_tile(aes(Pop1, Pop2, fill = jk_mean)) +
    facet_wrap(~T1_2)

#ggsave(file = "test.pdf", q, width = 10, height = 7)
#ggsave(file = "test.png", q, width = 10, height = 7)