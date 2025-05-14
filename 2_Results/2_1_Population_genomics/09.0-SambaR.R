
# Test whether working on HPC or laptop and set working directory accordingly
# Laptop test
if(grepl(getwd(), pattern = "C:/Users/mbzcp2/")){
  dir.path <-"C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/"
  plot.dir <- "C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/SambaR"
}
# HPC test
if(grepl(getwd(), pattern = "/gpfs01/home/mbzcp2")){
  dir.path <-"/gpfs01/home/mbzcp2/data/sticklebacks/"
  plot.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/SambaR"
}

## Set up SambaR
dir.create(paste0(plot.dir,"/Paired_with_outgroup"))
setwd(paste0(plot.dir,"/Paired_with_outgroup"))
source("/gpfs01/home/mbzcp2/R/x86_64-pc-linux-gnu-library/4.2/SambaR-master/SAMBAR_v1.10.txt")
## Install packages
getpackages()

## Read in data set 
importdata(inputprefix=paste0(dir.path,"vcfs/plink/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.recode"),
             samplefile=paste0(plot.dir, "/pop_file_Population.txt"))

filterdata(indmiss=0.5,snpmiss=0.2,min_mac=2 ,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=1000,nchroms=NULL,silent=TRUE) 

calcdiversity(nrsites=NULL,legend_cex=2.5)

findstructure(Kmax=6,add_legend=TRUE, legend_pos='bottomright', legend_cex=3, symbol_size=3) 

calckinship()
inferdemography()

calcdistance(nchroms=NULL) 
backupdata("Paired_with_outgroup")

## Clear out
rm(list = ls())

# Test whether working on HPC or laptop and set working directory accordingly
# Laptop test
if(grepl(getwd(), pattern = "C:/Users/mbzcp2/")){
  dir.path <-"C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/"
  plot.dir <- "C:/Users/mbzcp2/OneDrive - The University of Nottingham/Sticklebacks/Species Pairs M1/results/SambaR"
}
# HPC test
if(grepl(getwd(), pattern = "/gpfs01/home/mbzcp2")){
  dir.path <-"/gpfs01/home/mbzcp2/data/sticklebacks/"
  plot.dir <- "/gpfs01/home/mbzcp2/data/sticklebacks/results/SambaR"
}

# Move to specific directory
dir.creat(paste0(plot.dir,"/Paired_populaitons"))
setwd(paste0(plot.dir,"/Paired_populaitons"))
source("/gpfs01/home/mbzcp2/R/x86_64-pc-linux-gnu-library/4.2/SambaR-master/SAMBAR_v1.10.txt")
## Install packages
getpackages()

## Read in data set 
importdata(inputprefix=paste0(dir.path,"vcfs/plink/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.recode"),
             samplefile=paste0(plot.dir, "/pop_file_Population.txt"))

filterdata(indmiss=0.5,snpmiss=0.2,min_mac=2 ,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=1000,nchroms=NULL,silent=TRUE) 

subset_pop(include_pops=c("CLAC", "CLAM", "DUIM", "DUIN","LUIB", "LUIM","OBSE", "OBSM"))

calcdiversity(nrsites=NULL,legend_cex=2.5)

findstructure(Kmax=6,add_legend=TRUE, legend_pos='bottomright', legend_cex=3, symbol_size=3) 

calcdistance(nchroms=NULL) 

calckinship()
inferdemography()

backupdata("Paired_populaitons")


