args <- commandArgs(trailingOnly = TRUE)
## Read in bam stats
bam.stats <- read.table(args[1], header = T)

## Read  in sample data
samples <- read.table(args[2], header = T, sep =  ",")

# Merge two datasets
bam.stats <- merge(bam.stats, samples, by.x = "sample", by.y = "individual")
# Remove NAs
bam.stats <- bam.stats[!is.na(bam.stats$Population),]
# Replace stream designated populations with normal freshwater resident names
bam.stats$Population[bam.stats$Population=="OLST"] <- "OLAV"
bam.stats$Population[bam.stats$Population=="TOST"] <- "TORM"

## Select popualtions of interest
key.pops <- c("OBSE", "OBSM", "OLAM", "OLAV", "CLAM", "CLAC", "LUIB", "LUIM","DUIM", "DUIN")

# Subset data to those specific populations
bam.stats <- bam.stats[bam.stats$Population%in%key.pops,]

## Make mn_coverage numeric
bam.stats$mn_coverage <- as.numeric(gsub("X","", bam.stats$mn_coverage))

top.cov.samples <- list()
for(loc in key.pops){
  bam_tmp <- bam.stats[bam.stats$Population==loc,]
  top.cov.samples[[loc]] <- bam_tmp[(order(bam_tmp$mn_coverage, decreasing = T)[1]),]
  #print(paste("Top coverage samples for ", loc, " are:", paste(top.cov.samples[[loc]]$sample), "with ", paste(top.cov.outgroup[[loc]]$mn_coverage)))
}

top.cov.samples <- do.call("rbind", top.cov.samples)[,c("sample", "bam_file", "Waterbody","Population", "mn_coverage")]

## Write out table
write.table(top.cov.samples, file = args[3], row.names = F, quote = F, col.names = FALSE, sep = "\t")

