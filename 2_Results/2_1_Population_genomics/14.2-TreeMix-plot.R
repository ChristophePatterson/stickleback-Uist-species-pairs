source("/gpfs01/home/mbzcp2/apps/treemix-1.13/src/plotting_funcs.R")

## Get commnad arguments
args <- commandArgs(trailingOnly = TRUE)

# Get data
run_name <- args[1]

pdf(file = paste0(run_name,".pdf"), width = 10, height = 10)
plot_tree(cex=0.8,paste0(run_name))
dev.off()
