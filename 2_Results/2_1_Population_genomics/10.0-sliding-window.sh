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
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
conda activate genomics-general-p3.13

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/sliding-window

## Use Population (waterbody + ecotype)
# grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
#    awk -F ',' '{ print $1 " " $10}' > $wkdir/results/sliding-window/pop_file.txt

# Just use Ecotype
grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $13}' > $wkdir/results/sliding-window/pop_file.txt
# Pop file just for Y
grep -f $wkdir/vcfs/$vcf_ver/male_samples.txt $wkdir/results/sliding-window/pop_file.txt > $wkdir/results/sliding-window/pop_file_males.txt

# Run sliding window script
# For autosomes
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
   -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_auto.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/sliding-window/pop_file.txt -p anad -p resi

# For PAR
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz \
   -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_PAR.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/sliding-window/pop_file.txt -p anad -p resi

# For X
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz \
   -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_X.csv -f phased --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/sliding-window/pop_file.txt -p anad -p resi

# For Y
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.geno.gz \
   -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_Y.csv -f phased --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_Y.txt -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/sliding-window/pop_file_males.txt -p anad -p resi


# Merge all Fst calculations into a single file
cat $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_auto.csv > $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv
awk FNR!=1 $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_PAR.csv >> $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv
awk FNR!=1 $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_X.csv >> $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv
awk FNR!=1 $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi_Y.csv >> $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv

## Deactivate conda
conda deactivate
# Load R
module load R-uoneasy/4.2.1-foss-2022a
## Plot results
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi"
