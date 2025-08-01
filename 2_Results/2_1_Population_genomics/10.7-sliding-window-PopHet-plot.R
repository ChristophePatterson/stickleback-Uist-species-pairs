library(tidyverse)
library(patchwork)
library(tidyverse)
library(ggnewscale)

args <- commandArgs(trailingOnly = TRUE)

# set path
my_bins <- args[1]

## Set order pops should be plotted
allpops <- c("CLAC" ,"DUIN","LUIB","OBSE",  "CLAM", "DUIM", "LUIM", "OBSM")

# Read in data
hetPops <- read_csv(paste0(my_bins,".csv"))

## Pivot pi
hetPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("pi_", colnames(hetPops))])], 
             cols = colnames(hetPops)[grep("pi_", colnames(hetPops))], values_to = "pi", names_prefix = "pi_", names_to = "Population")
## Pivot dxy
dxyPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("dxy_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("dxy_", colnames(hetPops))], values_to = "dxy", names_prefix = "dxy_", names_to = "Population")
## Pivot fst
FstPops_long <- pivot_longer(hetPops[, c("scaffold","start","end","mid", "sites",colnames(hetPops)[grep("Fst_", colnames(hetPops))])], 
                             cols = colnames(hetPops)[grep("Fst_", colnames(hetPops))], values_to = "Fst", names_prefix = "Fst_", names_to = "Population")

## Summarise dxy
dxyPops_long_summary <- dxyPops_long %>%
  group_by(by = Population) %>%
  mutate(Pop1 = str_split_i(Population, "_", 1),
         Pop2 = str_split_i(Population, "_", 2)) %>%
  summarise(dxy.mn = mean(dxy, na.rm =T), dxy.sd = sd(dxy), 
            Pop1 = Pop1[1], Pop2 = Pop2[1])

## Order pop1 and pop2 so they are the same order as that in allpops
dxyPops_long_summary <- apply(dxyPops_long_summary, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops))]
}
      ) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(dxyPops_long_summary[,1:3])

## Summarise Fst
FstPops_long_summary <- FstPops_long %>%
  group_by(by = Population) %>%
  mutate(Pop1 = str_split_i(Population, "_", 1),
         Pop2 = str_split_i(Population, "_", 2)) %>%
  summarise(Fst.mn = mean(Fst, na.rm =T), Fst.sd = sd(Fst, na.rm = T), 
            Pop1 = Pop1[1], Pop2 = Pop2[1])

## Order pop1 and pop2 so they are the same order as that in allpops
FstPops_long_summary <- apply(FstPops_long_summary, MARGIN = 1, function(x) {
  c(x["Pop1"], x["Pop2"])[order(match(c(x["Pop1"], x["Pop2"]), allpops), decreasing = T)]
}
) %>%
  t() %>%
  as.data.frame() %>%
  cbind.data.frame(FstPops_long_summary[,1:3])

## Set factor level and reset POp1 and Pop2
dxyPops_long_summary$Pop1 <- factor(dxyPops_long_summary$V1,levels =  allpops)
dxyPops_long_summary$Pop2 <- factor(dxyPops_long_summary$V2,levels =  allpops)

## Set factor level and reset POp1 and Pop2
FstPops_long_summary$Pop1 <- factor(FstPops_long_summary$V1,levels =  allpops)
FstPops_long_summary$Pop2 <- factor(FstPops_long_summary$V2,levels =  allpops)

## plot
p <- ggplot() +
  geom_tile(data = dxyPops_long_summary, aes(Pop1, Pop2, fill = dxy.mn)) +
  scale_fill_gradient(low = "white", high = "deepskyblue",name = "Dxy") +
  new_scale_fill() +
  geom_tile(data = FstPops_long_summary, aes(Pop1, Pop2, fill = Fst.mn)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "white", high = "orange",name = "Fst") +
geom_vline(xintercept = 4.5) +
  geom_hline(yintercept = 4.5) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "grey")) +
  ggtitle("Dxy and Fst - Lagoon species pairs")


ggsave(paste0(my_bins, ".png"), p, height = 10, width = 10)
ggsave(paste0(my_bins, ".pdf"), p, height = 10, width = 10)

