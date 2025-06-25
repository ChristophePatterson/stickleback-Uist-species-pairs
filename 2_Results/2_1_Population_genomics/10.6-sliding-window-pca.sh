#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
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

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/
mkdir -p $outdir

########################
  # LEA - PCA & SNMF  # 
########################

## Run variables
vcf=$wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz
wndsize=500000
wndslid=250000

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
        $vcf $vcf_ver $wndsize $wndslid

