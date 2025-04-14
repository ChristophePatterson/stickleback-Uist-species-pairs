#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20g
#SBATCH --time=18:00:00
#SBATCH --job-name=sliding-window
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/sliding-window

## Use Population (waterbody + ecotype)
# grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
#    awk -F ',' '{ print $1 " " $10}' > $wkdir/results/sliding-window/pop_file.txt

# Just use Ecotype
grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
    awk -F ',' '{ print $1 " " $13}' > $wkdir/results/sliding-window/pop_file.txt


### Run sliding window script
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
    -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv -f phased -T 12 --popsFile $wkdir/results/sliding-window/pop_file.txt -p anad -p resi

cd $wkdir/results/sliding-window/
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/10.1-sliding-window-plot.R "sliding_window_w25kb_s5kb_m1_Panad_resi"

python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
    -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Pfw_anad.csv -f phased -T 12 --popsFile $wkdir/results/sliding-window/pop_file.txt -p fw -p anad

cd $wkdir/results/sliding-window/
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/10.1-sliding-window-plot.R "sliding_window_w25kb_s5kb_m1_Pfw_anad"


