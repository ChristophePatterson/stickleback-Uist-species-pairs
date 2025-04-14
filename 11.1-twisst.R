source("~/apps/twisst/plot_twisst.R")

## Get commnad arguments
args <- commandArgs(trailingOnly = TRUE)

# Get data
run_name <- args[1]

#weights file with a column for each topology
weights_file <- paste0(run_name, ".weights.tsv.gz")

# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
window_data_file <- paste0(run_name, ".data.tsv")

################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


pdf(file = paste0(run_name, ".summary.pdf"))
plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
dev.off()

pdf(file = paste0(run_name, ".all_typo.pdf"), width  = 10, height = 40)
plot.twisst(twisst_data, mode=1, show_topos=TRUE)
dev.off()




