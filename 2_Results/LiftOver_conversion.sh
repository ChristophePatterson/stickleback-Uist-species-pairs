## load liftover
conda activate LiftOver-env

##  Working directory
wkdir=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/)

# Unzip liftOver v1 to DUKE chain file
zcat $wkdir/gasAcu1.To.rabsTHREEspine.over.chain.gz > $wkdir/gasAcu1.To.rabsTHREEspine.over.chain.txt

# # # # # # # # # # # # # # # # # # #
# Convert gff from verion 5 to DUKE #
# # # # # # # # # # # # # # # # # # #

## Replace SeqID with Chromsome in R
# Unzip and convert to bed
zcat $wkdir/stickleback_v5_ensembl_genes.gff3.gz | awk '{print $1,$4,$5,$3,$9}' > $wkdir/stickleback_v5_ensembl_genes.bed

# Map v5 to v1
liftOver $wkdir/stickleback_v5_ensembl_genes.bed $wkdir/v5_to_v1_withChrUn.chain.txt \
    $wkdir/stickleback_v1_ensembl_genes.bed $wkdir/stickleback_v1_ensembl_genes_unlifted.bed -multiple -minMatch=0.01 bedPlus=3

# Tally up the number of unlifted files
grep '#' $wkdir/stickleback_v1_ensembl_genes_unlifted.bed | sort | uniq -c > $wkdir/stickleback_v1_ensembl_genes_unlifted_tally.txt
# Tally up the number of unlifted genes
grep -w "gene" -B 1 $wkdir/stickleback_v1_ensembl_genes_unlifted.bed | grep '#' | sort | uniq -c > $wkdir/stickleback_v1_ensembl_genes_unlifted_tally.txt

## v1 to DUKE
liftOver $wkdir/stickleback_v1_ensembl_genes.bed $wkdir/gasAcu1.To.rabsTHREEspine.over.chain.txt \
    $wkdir/stickleback_DUKE_ensembl_genes.bed $wkdir/stickleback_DUKE_ensembl_genes_unlifted.bed -multiple -minMatch=0.01  bedPlus=4

# Tally up the number of unlifted files
grep '#' $wkdir/stickleback_DUKE_ensembl_genes_unlifted.bed | sort | uniq -c > $wkdir/stickleback_DUKE_ensembl_genes_unlifted_tally.txt

# Count how many genes are not carryed over and why
grep -w "gene" -B 1 /gpfs01/home/mbzcp2/data/sticklebacks/genomes/stickleback_DUKE_ensembl_genes_unlifted.bed | grep '#' | sort | uniq -c

# # # # # # # # # # # # # #
### Jones et al to DUKE ###
# # # # # # # # # # # # # #
# Convert Jones et al CSS to DUKE genome coordinates

## change working directory
cd $wkdir/Prior_gasAcu-results

## Print bed and add unique name to each row
awk '{print $1,$2,$3,"j"NR}' Jones-et-al-2012-CSS-02-v4.bed > Jones-et-al-2012-CSS-02-v4_wID.bed

# v4 to v5
liftOver Jones-et-al-2012-CSS-02-v4_wID.bed ../v4_to_v5.chain.txt \
    Jones-et-al-2012-CSS-02-v5.bed Jones-et-al-2012-CSS-02-v5_unlifted.bed bedPlus=4 -multiple -minMatch=0

# v5 to v1
liftOver Jones-et-al-2012-CSS-02-v5.bed ../v5_to_v1.chain.txt \
    Jones-et-al-2012-CSS-02-v1.bed Jones-et-al-2012-CSS-02-v1_unlifted.bed bedPlus=4 -multiple -minMatch=0

# Change chr names
awk '{gsub(/group/, "chr");print}' Jones-et-al-2012-CSS-02-v1.bed > Jones-et-al-2012-CSS-02-v1_RN.bed

# v1 to DUKE
liftOver Jones-et-al-2012-CSS-02-v1_RN.bed ../gasAcu1.To.rabsTHREEspine.over.chain.txt \
    Jones-et-al-2012-CSS-02-vDUKE.bed Jones-et-al-2012-CSS-02-vDUKE_unlifted.bed bedPlus=4 -multiple -minMatch=0

# # # # # # # # # # # # # # # # # # # # # # # # # #
## Conversion of Roberts et al EcoPeaks to DUKE ###
# # # # # # # # # # # # # # # # # # # # # # # # # #
cd $wkdir/Prior_gasAcu-results

# v4 to v5
## Print bed and add unique name to each row
awk '{print $1,$2,$3,"r"NR}' Roberts-et-al-2021-Specific-EcoPeaks-v4.bed > Roberts-et-al-2021-Specific-EcoPeaks-v4_wID.bed
liftOver Roberts-et-al-2021-Specific-EcoPeaks-v4_wID.bed ../v4_to_v5.chain.txt \
    Roberts-et-al-2021-Specific-EcoPeaks-v5.bed Roberts-et-al-2021-Specific-EcoPeaks-v5_unlifted.bed bedPlus=4 -multiple minMatch=0

# v5 to v1
liftOver Roberts-et-al-2021-Specific-EcoPeaks-v5.bed ../v5_to_v1.chain.txt \
    Roberts-et-al-2021-Specific-EcoPeaks-v1.bed Roberts-et-al-2021-Specific-EcoPeaks-v1_unlifted.bed bedPlus=4 -multiple minMatch=0

# Change chr names
awk '{gsub(/group/, "chr");print}' Roberts-et-al-2021-Specific-EcoPeaks-v1.bed > Roberts-et-al-2021-Specific-EcoPeaks-v1_RN.bed

# v1 to DUKE
liftOver Roberts-et-al-2021-Specific-EcoPeaks-v1_RN.bed ../gasAcu1.To.rabsTHREEspine.over.chain.txt \
    Roberts-et-al-2021-Specific-EcoPeaks-vDUKE.bed Roberts-et-al-2021-Specific-EcoPeaks-vDUKE_unlifted.bed bedPlus=4 -multiple minMatch=0

## deactivate liftOver
conda deactivate

### Plot which regions were not lifted over
# Load R
module load R-uoneasy/4.2.1-foss-2022a

cd $wkdir

# Run R script to visualise conversions
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/LiftOver_conversion.R
