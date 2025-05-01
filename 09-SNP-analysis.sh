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
#SBATCH --job-name=SNP-PCA-LEA
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/09-SNP-analysis.R

conda activate genomics-general-p3.13

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/popgen
mkdir -p $outdir

## Create popfile
grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $10}' > $outdir/pop_file_Population.txt

python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis indHet -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
   -o $outdir/sliding_window_w25kb_s5kb_m1_Popgen.csv -f phased -T $SLURM_CPUS_PER_TASK \
   --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 

python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis hapStats -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
   -o $outdir/sliding_window_w25kb_s5kb_m1_hapStats.csv -f phased -T $SLURM_CPUS_PER_TASK \
   --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 


conda deactivate