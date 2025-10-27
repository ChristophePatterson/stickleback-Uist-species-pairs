#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH --time=18:00:00
#SBATCH --job-name=SNP-PCA-LEA
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
species=stickleback
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver
mkdir -p $outdir

mkdir -p $wkdir/vcfs/$vcf_ver/plink

## ## Plink/Sambar cant have samples with "_" - replace with "-"
## bcftools query -l $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz > $wkdir/vcfs/$vcf_ver/plink/samples.txt
## sed s/_/-/ $wkdir/vcfs/$vcf_ver/plink/samples.txt > $wkdir/vcfs/$vcf_ver/plink/samples_recode.txt
## # Reheader samples
## bcftools reheader -s $wkdir/vcfs/$vcf_ver/plink/samples_recode.txt $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz > $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.vcf.gz
## # deactivate bcftools
## conda deactivate
## 
## module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig
## ## Convert to Plink format to inlcude input into SambaR
## plink -vcf  $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader
## plink --file $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.recode
## 
## ## Create input popfile for SambaR
## echo -e "name\tpop" > $wkdir/vcfs/$vcf_ver/plink/pop_file_Population.txt
## grep -f $wkdir/vcfs/$vcf_ver/plink/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
##     awk -F ',' -v OFS='\t' '{ print $1, $10}' | sed s/NA/Lubec/ | sed s/_/-/ >> $wkdir/vcfs/$vcf_ver/plink/pop_file_Population.txt

# Deactivate conda
module purge

module load R-uoneasy/4.2.1-foss-2022a
# ## Run SambaR
# Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09.0-SambaR.R \
#   $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.recode

################################
  # Genomic sex determination  # 
################################

##### Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/Sex_determination.R

########################
  # LEA - PCA & SNMF  # 
########################

## Custom analysis

## Masked data
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09-SNP-analysis.R \
       $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.geno $vcf_ver
#### 
#### # Unmasked data
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09-SNP-analysis.R \
       $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000.geno $vcf_ver


## Conduct analysis on hard filtered FST outliers
module purge
conda activate bcftools-env

# Set FST outlier threshold
FstUpper="0.1"

## Get regions of genome that are FST outliers in bed format
awk -F "," -v OFS='\t' -v Fst=$FstUpper 'NR!=1 && $9!="nan" && $9 >= (Fst + 0) { print $1, $2, $3, $9}' \
    $wkdir/results/$vcf_ver/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_auto.csv > \
    $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.FST.$FstUpper.outliers.bed


## Filter out snps from these regions
bcftools view -T ^$wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.FST.$FstUpper.outliers.bed \
    -O z -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.noFST$FstUpper.vcf.gz \
    $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Filter to get 1 snp per 1kb
bcftools +prune -n 1 -N rand -w 1000bp -O v -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.noFST$FstUpper.rand1000.vcf.gz \
    $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.noFST$FstUpper.vcf.gz

# Load R module
module purge
module load R-uoneasy/4.2.1-foss-2022a

# Convert to geno format
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/08.1-vcf2geno.R \
    $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.noFST$FstUpper.rand1000.vcf.gz \
    $vcf_ver

# Run LEA analysis
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09-SNP-analysis.R \
        $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.noFST$FstUpper.rand1000.geno $vcf_ver

