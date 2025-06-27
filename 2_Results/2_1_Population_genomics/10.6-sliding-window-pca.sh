#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=18:00:00
#SBATCH --job-name=sliding-window-pca
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware_HWEPops_MQ10_BQ20

########################
  # LEA - PCA & SNMF  # 
########################

## Run variables
vcf=${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.r1000
wndsize=1000000
wndslid=500000

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/$vcf
mkdir -p $outdir


Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
        $wkdir/vcfs/$vcf_ver/$vcf.vcf.gz $vcf_ver $wndsize $wndslid
        
## Analysis where population can be changed (specifically removing CLAC)
run_name=(noCLAC)
outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/$vcf.$run_name
mkdir -p $outdir

# Get names of samples in vcf
bcftools query -l $wkdir/vcfs/$vcf_ver/$vcf.vcf.gz > $outdir/samples_in_vcf.txt
# Create file of waterbodies that are of interes
echo -e "OBSE\nDUIN\nLUIB" > $outdir/Pops_interest.txt
# Get sample information and subset to those water bodies
awk -F ',' -v OFS='\t' '{ print $1, $13 ,$9}' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
         grep -f $outdir/samples_in_vcf.txt | \
         grep -f $outdir/Pops_interest.txt > $outdir/pop_file.txt
# Get subset list of sample files
awk '{print $1}' $outdir/pop_file.txt > $outdir/samples.txt

# Remove CLAC samples and remove nolonger variable snps
bcftools view -S $outdir/samples.txt --min-ac 2:minor -O z -o $outdir/$vcf.$run_name.vcf.gz $wkdir/vcfs/$vcf_ver/$vcf.vcf.gz

# Run analysis again
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
        $outdir/$vcf.$run_name.vcf.gz $vcf_ver $wndsize $wndslid

rm $outdir/$vcf.$run_name.vcf.gz






