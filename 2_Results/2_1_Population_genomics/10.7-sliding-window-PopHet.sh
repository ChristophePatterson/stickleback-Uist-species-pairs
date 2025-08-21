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
#SBATCH --job-name=sliding-window-PopHet
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
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
output_dir=($wkdir/results/$vcf_ver/sliding-window/indPops)
mkdir -p ${output_dir}

## Use Population (waterbody + ecotype)
grep -w -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
   awk -F ',' '{ print $1 " " $10}' > ${output_dir}/pop_file.txt

# Just use Ecotype
# grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
#    awk -F ',' -v OFS='\t' '{ print $1, $13}' > ${output_dir}/pop_file.txt

# Female file just for X
grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt ${output_dir}/pop_file.txt > ${output_dir}/pop_file_females.txt
# Pop file just for Y
grep -w -f $wkdir/vcfs/$vcf_ver/male_samples.txt ${output_dir}/pop_file.txt > ${output_dir}/pop_file_males.txt

# Run PopDist window script
# For autosomes
python ~/apps/genomics_general/popgenWindows.py -w 100000 -s 100000 -m 1 --analysis popDist popPairDist indHet -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
   -o ${output_dir}/sliding_window_w100kb_s100kb_m1_PopPair_auto.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile ${output_dir}/pop_file.txt -p DUIN -p OBSE -p LUIB -p CLAC -p DUIM -p OBSM -p LUIM -p CLAM

## Deactivate conda
conda deactivate
# Load R
module load R-uoneasy/4.2.1-foss-2022a
## Plot results
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.7-sliding-window-PopHet-plot.R "${output_dir}/sliding_window_w100kb_s100kb_m1_PopPair_auto"
