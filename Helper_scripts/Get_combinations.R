## Code for getting all unique combination
## args <A text file containing values on each row that you want each unique combination off>
args <- commandArgs(trailingOnly=T)
values <- readLines(args[1])
ncomb <- as.numeric(args[2])
comb.values <- t(combn(unique(values), ncomb))
write.table(file = paste0(args[1],"_combn.txt"), comb.values, row.names = F, quote = F, col.names = F, sep = "\t")
