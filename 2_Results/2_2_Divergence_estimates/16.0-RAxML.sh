#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=6g
#SBATCH --time=24:00:00
#SBATCH --job-name=RaxML
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

## Load RAxML
module load raxml-ng-uoneasy/1.2.0-GCC-12.3.0

# Create output directory
output_dir=$wkdir/results/$vcf_ver/phylos/RaxML
mkdir -p $output_dir

## Run RAxML
SNP_library=stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000
phy_file="$SNP_library.phy"

# Change into output directory
cd $output_dir

rm *${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS*

# Run RAxML
pwd
echo $phy_file

# raxml-ng --check --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --threads 30 --prefix GTGTR4_G_ASC_LEWIS
raxml-ng --all --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --threads 30 --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS
