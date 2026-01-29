library(tidyverse)
library(patchwork)
library(ggnewscale)

# Colourblind palette
cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")

# Model name
model_dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/lochs_mono_array_folded_nCDS_nHFst_r10000/models_lochs_mono_array_folded_nCDS_nHFst_r10000_all_results"
setwd(model_dir)
model_name <- basename(model_dir)

# Get file names once
files <- list.files(pattern = "\\.bestlhoods")

# Derive run names
runs <- sub("\\.bestlhoods", "", files)

# Read all files once
runs.list <- lapply(files, read.table, header = TRUE)

# Get union of all column names
var.names <- Reduce(union, lapply(runs.list, colnames))

# Pre-allocate result matrix
runs.data <- matrix(
  NA,
  nrow = length(runs.list),
  ncol = length(var.names),
  dimnames = list(runs, var.names)
)

i <- 1
# Fill matrix
for (i in seq_along(runs.list)) {
  tmp <- runs.list[[i]]
  idx <- match(colnames(tmp), var.names)
  runs.data[i, idx] <- unlist(c(tmp[1, ]))
}
runs.data <- as.data.frame(runs.data)

runs.data$run <- row.names(runs.data)

# Tidy up data frame to contain metadata for each run
full.model.data <- as.data.frame(runs.data) %>%
  mutate(pop0 = str_split_i(run, "-", 1),
         pop1 = str_split_i(run, "-", 2),
         model = str_split_i(run, "-", 3),
         foldtype = str_split_i(run, "-", 4),
         SFS = str_split_i(run, "-", 5))

## HOw many vars are in each model
full.model.data$N <- apply(full.model.data, MARGIN = 1, function(x) sum(!is.na(x[var.names[!grepl("Lhood", var.names)]])))

full.model.data <- full.model.data %>%
  mutate(AIC = (2*N)-(2*MaxEstLhood))

# Get population names
pop.order <- c("CLAC","CLAM", "DUIN", "DUIM", "LUIB","LUIM", "OBSE", "OBSM")
# Factor populations
full.model.data$pop0 <- factor(full.model.data$pop0,levels = pop.order)
full.model.data$pop1 <- factor(full.model.data$pop1,levels = pop.order)
# Factor models by complexity
model.complex.order <- full.model.data$model[!duplicated(full.model.data$model)][order(full.model.data$N[!duplicated(full.model.data$model)])]
full.model.data$model <-  factor(full.model.data$model, levels = model.complex.order)

# Pivot to long format
full.model.data.long <- pivot_longer(full.model.data, cols = c("AIC",var.names)) %>%
  mutate(is.div = grepl("DIV|TMIG|TNE", name),
         is.popNe = !grepl("DIV|TMIG|Lhood|AIC|MIG|TNE",name),
         is.NPOP = grepl("NPOP",name),
         is.MIG = grepl("MIG",name)&!grepl("TMIG",name))

ggplot(full.model.data.long[full.model.data.long$is.popNe,]) +
  geom_boxplot(aes(name, value/2, col = model)) +
  facet_wrap(~pop0*model, scales = "free_x", drop = T) +
  scale_y_continuous(labels = function(x) x/1e6, name = "Ne (Millions)")

ggplot(full.model.data.long[full.model.data.long$is.NPOP,]) +
  geom_boxplot(aes(name, value/2, col = model)) +
  geom_jitter(aes(name, value/2, col = model), height = 0) +
  facet_wrap(~pop0*model, scales = "free_x", drop = T) +
  scale_y_continuous(labels = function(x) x/1e6, name = "Ne (Millions)")

ggplot(full.model.data.long[full.model.data.long$is.div,]) +
  geom_boxplot(aes(name, value, col = model)) +
  facet_wrap(~pop0) +
  # scale_y_log10(labels = function(x) x/1e6, name = "Years (Millions)")
  scale_y_continuous(labels = function(x) x/1e6, name = "Millions of Generations (1 gen ~ year)")

ggplot(full.model.data.long[full.model.data.long$is.MIG,]) +
  geom_boxplot(aes(name, value, col = model)) +
  facet_wrap(~pop0) +
  scale_y_log10()

ggplot(full.model.data) +
  geom_boxplot(aes(model, AIC, col = model), outliers = F) +
  geom_jitter(aes(model, AIC, col = model), height = 0) +
  facet_wrap(~pop0, scales = "free_y")

ggplot(full.model.data) +
  geom_jitter(aes(model, MaxEstLhood-MaxObsLhood, col = model), height = 0) +
  facet_wrap(~pop0, scales = "free_y") + 
  scale_color_manual(values = cbPalette, name = "Loch") +
  theme_bw()

table(full.model.data$pop0, full.model.data$SFS)

# Calculate delta AIC
delta.AIC <- full.model.data %>%
  group_by(pop0) %>%
  group_by(SFS, .add = T) %>%
  reframe(min.AIC = min(AIC),
            min.model = model[which(AIC==min(AIC))],
            delta.AIC.Iso = min.AIC-AIC[model=="Iso"],
            delta.AIC.IsoMig = min.AIC-AIC[model=="IsoMig"],
            delta.AIC.IsoSC = min.AIC-AIC[model=="IsoSC"],
            delta.AIC.IsoNeC = min.AIC-AIC[model=="IsoNeC"]) %>%
  pivot_longer(cols = c("delta.AIC.Iso", "delta.AIC.IsoMig", "delta.AIC.IsoSC", "delta.AIC.IsoNeC"),
               names_to = "model",values_to = "d.AIC") %>%
  mutate(model = factor(gsub("delta.AIC.","",model), levels = model.complex.order))

table(delta.AIC$pop0, delta.AIC$min.model)/4

# Summarise best model counts
best.model.df <- delta.AIC %>%
  group_by(pop0, min.model) %>%
  mutate(min.model = factor(min.model, levels = model.complex.order)) %>%
  summarise(best.model.count = n()/4, ) %>%
  rename(model = min.model)
# Plot delta AIC
p.AIC <- ggplot(delta.AIC, aes(col = pop0)) +
  geom_boxplot(aes(model, d.AIC), width = 0.5, outliers = F) +
  geom_jitter(aes(model, d.AIC, fill = (d.AIC == 0)),
              size = 2, shape = 21, height = 0, width = 0.2) +
  geom_text(data = best.model.df, aes(model, 0, label = best.model.count), col = "black", vjust = -1) +
  facet_wrap(~pop0, scales="free_y") +
  scale_y_continuous(name = "ΔAIC", limits = c(-25, 2),, oob = scales::squish) +
  scale_color_manual(name = "Loch", values = cbPalette) +
  scale_fill_manual(values = c("white","black"), name = "Min AIC") +
  theme_bw()

p.AIC

p <- ggplot(full.model.data, aes(col = pop0)) +
  geom_violin(aes(x = TDIV, y = as.numeric(pop0)+0.5), width = 0.2) +
  geom_segment(aes(x = TDIV, y = pop0, yend = pop1), linewidth = 1, alpha = 0.25) +
  geom_segment(aes(x = TDIV, xend = 0, y = pop0, yend = pop0, linewidth = NPOP0/2), alpha = 0.25) +
  geom_segment(aes(x = TDIV, xend = 0, y = pop1, yend = pop1, linewidth = NPOP1/2), alpha = 0.25) +
  geom_segment(aes(x = TMIG, y = pop0, yend = pop1), linewidth = 1, alpha = 0.25) +
  facet_wrap(~model,ncol = 3, scales = "free_x") +
  scale_y_discrete(drop = F, name = "Population") +
  scale_x_reverse(expand = c(0.15, 0), name = "Generations (~1 year)") +
  scale_linewidth_continuous(labels = function(x) x/1e6, name = "Ne (Millions)") +
  scale_color_discrete(name = "Loch") +
  theme_bw()

p

# Create summary stats
median.model <- full.model.data %>%
  group_by(pop0, model, .add = T) %>%
  summarise(pop1 = first(pop1),
            md.TDIV = median(TDIV, na.rm = T),
            md.TMIG = median(TMIG, na.rm = T),
            md.MIG01 = median(MIG01, na.rm = T),
            md.MIG10 = median(MIG10, na.rm = T),
            md.NPOP0 = median(NPOP0, na.rm = T),
            md.NPOP1 = median(NPOP1, na.rm = T),
            md.ANC0 = median(ANC0, na.rm = F),
            md.ANC1 = median(ANC1, na.rm = F),
            md.POP0TNE = median(POP0TNE, na.rm = F),
            md.POP1TNE = median(POP1TNE, na.rm = F),
            sd.TDIV = sd(TDIV, na.rm = T))

# Manipulate NPOP and ANC for plotting change in Ne
median.model$md.POP0TNE[is.na(median.model$md.POP0TNE)] <- 0
median.model$md.POP1TNE[is.na(median.model$md.POP1TNE)] <- 0
median.model$md.ANC0[is.na(median.model$md.ANC0)] <- median.model$md.NPOP0[is.na(median.model$md.ANC0)]
median.model$md.ANC1[is.na(median.model$md.ANC1)] <- median.model$md.NPOP1[is.na(median.model$md.ANC1)]

# Plot with medians
p <- ggplot(full.model.data, aes(col = pop0)) +
  geom_violin(aes(x = TDIV, y = as.numeric(pop0)+0.5), width = 0.5) +
  geom_violin(aes(x = TMIG, y = pop0), width = 0.5) +
  geom_violin(aes(x = TMIG, y = pop1), width = 0.5) +
  geom_point(data = median.model, aes(x = md.TDIV, y = as.numeric(pop0)+0.5)) +
  geom_point(data = median.model, aes(x = md.TMIG, y = pop0)) +
  geom_point(data = median.model, aes(x = md.TMIG, y = pop1)) +
  geom_segment(data = median.model, aes(x = md.TDIV, y = pop0, yend = pop1), linewidth = 1, alpha = 0.5) +
  geom_segment(data = median.model, aes(x = md.TDIV, xend = md.POP0TNE, y = pop0, yend = pop0, linewidth = md.ANC0/2), alpha = 0.5) +
  geom_segment(data = median.model, aes(x = md.TDIV, xend = md.POP1TNE, y = pop1, yend = pop1, linewidth = md.ANC1/2), alpha = 0.5) +
  geom_segment(data = median.model, aes(x = md.POP0TNE, xend = 0, y = pop0, yend = pop0, linewidth = md.NPOP0/2), alpha = 0.5) +
  geom_segment(data = median.model, aes(x = md.POP1TNE, xend = 0, y = pop1, yend = pop1, linewidth = md.NPOP1/2), alpha = 0.5) +
  #geom_segment(data = median.model, aes(x = md.TDIV/2, y = pop0, yend = pop1), linewidth = 1, linetype = "dashed") +
  #geom_segment(data = median.model, aes(x = md.TDIV/2, y = pop1, yend = pop0), linewidth = 1, linetype = "dashed") +
  geom_segment(data = median.model, aes(x = md.TMIG, y = pop0, yend = pop1), linewidth = 1, arrow = arrow(length = unit(10,units = "points"))) +
  geom_segment(data = median.model, aes(x = md.TMIG, y = pop1, yend = pop0), linewidth = 1, arrow = arrow(length = unit(10,units = "points"))) +
  geom_text(data = median.model, aes(x = md.TMIG, y = pop0, label = round(md.NPOP0*md.MIG01)), vjust = 1) +
  geom_text(data = median.model, aes(x = md.TMIG, y = pop1, label = round(md.NPOP1*md.MIG10)), vjust = -1) +
  geom_text(data = median.model, aes(x = 0, y = pop0, label = round(md.NPOP0/2)), hjust = -0.1) +
  geom_text(data = median.model, aes(x = 0, y = pop1, label = round(md.NPOP1/2)), hjust = -0.1) +
  facet_wrap(~model,ncol = 2, scales = "free_x") +
  scale_y_discrete(drop = F, name = "Population") +
  scale_x_reverse(expand = c(0.23, 0), labels = function(x) x, name = "Generations (~1 year)") +
  scale_linewidth_continuous(labels = function(x) x/1e6, name = "Ne (Millions)") +
  scale_color_manual(name = "Loch", values = cbPalette) +
  theme_bw()
p

ggsave(paste0(model_name,".png"), p / p.AIC + plot_layout(heights = c(2,1)), width = 10, height = 15)
