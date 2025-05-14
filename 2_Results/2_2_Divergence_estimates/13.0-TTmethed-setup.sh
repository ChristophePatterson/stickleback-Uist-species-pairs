#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=48:00:00
#SBATCH --job-name=twisst
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/twisst-%x-%j.out


## Load modules
module load R-uoneasy/4.2.1-foss-2022a

## Set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback

## Create directories
mkdir -p $wkdir/results/TTmethod
## Create folder for sample vcfs
mkdir -p $wkdir/results/TTmethod/vcfs

##  Create data set of which samples to use
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/13.1-TTmethod-setup.R \
    /gpfs01/home/mbzcp2/data/sticklebacks/bams/bamstats/QC/clean_bams/Multi-Bam-QC/global_raw_report_custom.txt \
    /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv \
    $wkdir/results/TTmethod/vcfs/top_coverage_samples.txt

## Run calling on all top covarage samples
sbatch /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/13.2-calling-top-cov.sh






