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
module load R-uoneasy/4.2.1-foss-2022a
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/SambaR
mkdir -p $outdir

## Plink/Sambar cant have samples with "_" - replace with "-"
bcftools query -l $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz > $outdir/samples.txt
sed s/_/-/ $outdir/samples.txt > $outdir/samples_recode.txt
# Reheader samples
bcftools reheader -s $outdir/samples_recode.txt $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.vcf.gz

## Convert to Plink format to inlcude input into SambaR
mkdir -p $wkdir/vcfs/plink
plink -vcf  $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader
plink --file $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.recode

## Create input popfile for SambaR
echo -e "name\tpop" > $outdir/pop_file_Population.txt
grep -f $outdir/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $10}' | sed s/NA/Lubec/ | sed s/_/-/ >> $outdir/pop_file_Population.txt

## Run SambaR
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/09.0-SambaR.R

## Custom analysis
# Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/09-SNP-analysis.R

conda activate genomics-general-p3.13

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/popgen
mkdir -p $outdir

## Create popfile
## grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
##     awk -F ',' -v OFS='\t' '{ print $1, $10}' > $outdir/pop_file_Population.txt
## 
## python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis indHet -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
##    -o $outdir/sliding_window_w25kb_s5kb_m1_Popgen.csv -f phased -T $SLURM_CPUS_PER_TASK \
##    --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 
## 
## python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis hapStats -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
##    -o $outdir/sliding_window_w25kb_s5kb_m1_hapStats.csv -f phased -T $SLURM_CPUS_PER_TASK \
##    --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 


conda deactivate