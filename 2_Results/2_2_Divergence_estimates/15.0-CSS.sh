#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH --time=18:00:00
#SBATCH --job-name=CSS
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
conda activate bcftools-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
output_dir=$wkdir/results/sliding-window/CSS
# Create 
mkdir -p $output_dir

## Script output to location of vcf and needs cleaner file name
# Define which vcf to use
vcf=/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/ploidy_aware/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz
## Copy over to outpute folder
cp $vcf $output_dir/stickleback.vcf.gz

# Query names of samples in VCF
bcftools query -l $vcf > $output_dir/samples.txt

## Create pop file for samples
grep -f $output_dir/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $13}' > $output_dir/pop_file.txt

## Deactivate conda
conda deactivate
# Load R
module load R-uoneasy/4.2.1-foss-2022a
## Plot results

## Remove prior results if they exist
if [ -f $output_dir/stickleback.gds ]; then
   echo "gds file already exists so removing"
   rm $output_dir/stickleback.gds 
else
   echo "gds file does not exists, running analysis."
fi


# Change to output directy
cd $output_dir
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/CSSm.R \
         stickleback.vcf.gz pop_file.txt 2500 500 10 locus mds "0.05"
