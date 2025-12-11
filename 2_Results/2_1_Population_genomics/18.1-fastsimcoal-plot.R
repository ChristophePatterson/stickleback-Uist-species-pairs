library(tidyverse)
library(patchwork)

## only currently works on uni laptop (not Ada)

# Set working directory
setwd("C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/Demographic/fsc_run/")

# Get list of pop pairs
pop.pairs <- data.frame(pops = gsub(".par", "", list.files()[grep(".par", list.files())]))

# Get population names
pop.pairs$pop1 <- str_split_i(pop.pairs$pops, "-", 1)
pop.pairs$pop2 <- str_split_i(pop.pairs$pops, "-", 2)

## Get best module for each pair
sim.limit <- do.call("rbind", apply(pop.pairs, MARGIN = 1, function(x) read.table(paste0(x["pops"], ".est"), header = F, skip = 5, nrows = 5,
                                                                                  col.names = c("int", "param", "shape", "min", "max", "print"))))

## Get best module for each pair
pop.pairs <- cbind(pop.pairs, 
                   do.call("rbind", apply(pop.pairs, MARGIN = 1, function(x) read.table(paste0(x["pops"], "/", x["pops"], ".bestlhoods"), header = T))))

## Dived N pop by 2 to convert haploid to diploid
pop.pairs$ANCSIZE <- pop.pairs$ANCSIZE/2
pop.pairs$NPOP1 <- pop.pairs$NPOP1/2
pop.pairs$NPOP2 <- pop.pairs$NPOP1/2

## Get all param calculations
besthoods.df <-  do.call("rbind", apply(pop.pairs, MARGIN = 1, 
                                        function(x) {
                                          besthood.sim.length <- length(readLines(paste0(x["pops"], "/", x["pops"], ".brent_lhoods")))-4 # Get number of simulations
                                          return(
                                            cbind.data.frame(INT = 1:besthood.sim.length, pops = x["pops"], pop1 = x["pop1"], pop2 = x["pop2"],read.table(paste0(x["pops"], "/", x["pops"], ".brent_lhoods"),
                                                                                                                                                          col.names = c("PARAM", "ANCSIZE","NPOP1","NPOP2","TDIV","MIG21","MIG12","MaxEstLhood"),
                                                                                                                                                          skip = 1, nrows = besthood.sim.length)
                                            ))
                                        }))
## Dived N pop by 2 to convert haploid to diploid
besthoods.df$ANCSIZE <- besthoods.df$ANCSIZE/2
besthoods.df$NPOP1 <- besthoods.df$NPOP1/2
besthoods.df$NPOP2 <- besthoods.df$NPOP1/2


## Plot of best pop sizes
ggplot(pop.pairs) +
  geom_hline(yintercept = sim.limit$min[sim.limit$param=="NPOP1"]/2) +
  geom_hline(yintercept = sim.limit$max[sim.limit$param=="NPOP1"]/2) +
  geom_point(aes(pop1, NPOP1)) +
  geom_point(aes(pop2, NPOP2))

## Plot of best migrates
ggplot(pop.pairs) +
  geom_hline(yintercept = sim.limit$min[sim.limit$param=="N1M21"]/sim.limit$min[sim.limit$param=="NPOP1"]) +
  geom_hline(yintercept = sim.limit$max[sim.limit$param=="N2M12"]/sim.limit$max[sim.limit$param=="NPOP2"]) +
  geom_point(aes(pop1, MIG12)) +
  geom_point(aes(pop2, MIG21))


ggplot(besthoods.df) +
  geom_hline(yintercept = sim.limit$min[sim.limit$param=="NPOP1"]/2) +
  geom_hline(yintercept = sim.limit$max[sim.limit$param=="NPOP1"]/2) +
  geom_point(aes(INT, NPOP1, col = pop1)) +
  geom_point(aes(INT, NPOP2, col = pop2)) +
  geom_point(data = pop.pairs, aes(max(besthoods.df$INT)*1.1, NPOP1, fill = pop1), size = 3, shape = 23) +
  geom_point(data = pop.pairs, aes(max(besthoods.df$INT)*1.1, NPOP2, fill = pop2), size = 3, shape = 23) +
  geom_text(data = pop.pairs, aes(max(besthoods.df$INT)*1.1, NPOP1, col = pop1, label = NPOP1), hjust =-0.2, show.legend = F) +
  geom_text(data = pop.pairs, aes(max(besthoods.df$INT)*1.1, NPOP2, col = pop2, label = NPOP2), hjust =-0.2, show.legend = F) +
  scale_x_continuous(limits = c(0, max(besthoods.df$INT)*1.2))

