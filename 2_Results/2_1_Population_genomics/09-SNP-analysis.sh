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

conda activate bcftools-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware_HWEPops_MQ10_BQ20

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/SambaR
mkdir -p $outdir

mkdir -p $wkdir/vcfs/$vcf_ver/plink

## Plink/Sambar cant have samples with "_" - replace with "-"
bcftools query -l $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz > $wkdir/vcfs/$vcf_ver/plink/samples.txt
sed s/_/-/ $wkdir/vcfs/$vcf_ver/plink/samples.txt > $wkdir/vcfs/$vcf_ver/plink/samples_recode.txt
# Reheader samples
bcftools reheader -s $wkdir/vcfs/$vcf_ver/plink/samples_recode.txt $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz > $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.vcf.gz
# deactivate bcftools
conda deactivate

module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig
## Convert to Plink format to inlcude input into SambaR
plink -vcf  $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader
plink --file $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/$vcf_ver/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.reheader.recode

## Create input popfile for SambaR
echo -e "name\tpop" > $wkdir/vcfs/$vcf_ver/plink/pop_file_Population.txt
grep -f $wkdir/vcfs/$vcf_ver/plink/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $10}' | sed s/NA/Lubec/ | sed s/_/-/ >> $wkdir/vcfs/$vcf_ver/plink/pop_file_Population.txt

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
        $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.vcf.gz $vcf_ver

# Unmasked data
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09-SNP-analysis.R \
        $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz $vcf_ver

#####################
  # fastStructure  # 
#####################

## conda create -n faststucture-p2 python=2
###### module purge
###### conda activate faststucture-p2
###### module load bcftools-uoneasy/1.18-GCC-13.2.0
###### module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig
###### 
###### outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/faststructure
###### mkdir -p $outdir
###### 
###### ## Create input file for fasstructure
###### inputfile=(stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000)
###### 
###### ## Remove non paired species 
###### bcftools view -S $wkdir/vcfs/${species}_subset_samples.txt $wkdir/vcfs/${inputfile}.vcf.gz | \
######     bcftools view --min-ac 2[minor] -O z -o $wkdir/vcfs/${inputfile}.rand1000.SpPair.vcf.gz
###### 
###### ## Convert to plink format
###### plink -vcf  $wkdir/vcfs/${inputfile}.rand1000.SpPair.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair
###### plink --file $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair.recode
###### 
###### ## Create popfile for plotting results
###### grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
######     awk -F ',' -v OFS='\t' '{ print $10 }' > $outdir/popfile.txt
###### 
###### ## Loop through running faststructure with K= 1 to 5
###### # Set i as 0
###### i=0
###### ## Start looop
###### while [ $i -ne 20 ]
###### do
######         ## Add 1 to K
######         i=$(($i+1))
######         echo "Running faststructre for K-$i"
######         ## Run structure
######         structure.py --input=$wkdir/vcfs/plink/${inputfile}.rand1000.SpPair.recode \
######             --format bed -K $i --output=$outdir/${species}_SpPair \
######             --cv=20 --full
###### 
######         ## Plot stucture make sure you have added code to distruct.py
######         # beforethe to the code "import matplotlib.pyplot as plot"
######         # #import matplotlib as mpl
######         # #mpl.use('Agg')
###### 
######         distruct.py --input=$outdir/${species}_SpPair \
######              -K $i --output=$outdir/${species}_SpPair_K${i}.pdf --popfile=$outdir/popfile.txt \
######              --title=${species}_SpPair_K${i}
###### 
###### done
###### 
###### conda deactivate
###### 
