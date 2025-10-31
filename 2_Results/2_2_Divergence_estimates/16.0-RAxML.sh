#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=93
#SBATCH --mem=25g
#SBATCH --time=167:00:00
#SBATCH --job-name=RAxML
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

## Load RAxML
module load raxml-ng-uoneasy/1.2.0-GCC-12.3.0

# Create output directory
output_dir=$wkdir/results/$vcf_ver/phylos/RAxML
mkdir -p $output_dir

## Run RAxML
SNP_library=stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000
phy_file="$SNP_library.phy"

# Change into output directory
cd $output_dir

# Run RAxML
pwd
echo $phy_file

# raxml-ng --check --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --threads $SLURM_CPUS_PER_TASK --prefix GTGTR4_G_ASC_LEWIS
# raxml-ng --parse --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --prefix GTGTR4_G_ASC_LEWIS
raxml-ng --all --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --bs-trees 500 --tree pars{10},rand{10} --threads $SLURM_CPUS_PER_TASK --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10

## Map bootstrap values on the ML tree
raxml-ng --support --tree ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bestTree --bs-trees ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bootstraps --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml --threads $SLURM_CPUS_PER_TASK

## Test convergence of bootstrap trees
raxml-ng --bsconverge --bs-trees ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bootstraps --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bsconverge --threads $SLURM_CPUS_PER_TASK --bs-cutoff 0.01


## module purge
## Load R
module load R-uoneasy/4.2.1-foss-2022a

# Run R plotting script
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/16.1-RAxML_plot.R \
    ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml
