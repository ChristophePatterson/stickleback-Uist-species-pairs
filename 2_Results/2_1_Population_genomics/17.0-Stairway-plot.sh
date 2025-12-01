#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=18:00:00
#SBATCH --array=1-8
#SBATCH --job-name=stairway_plot
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc

# Load easySFS
conda activate easySFS-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

## Output
output_dir=($wkdir/results/$vcf_ver/Stairway)
mkdir -p $output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
# echo -e "CLAC\nCLAM\nOBSE\nOBSM\nDUIN\nDUIM\nLUIB\nLUIM" > ${output_dir}/pop_list.txt
## Get population that equals slurm array
pop=$(awk -v slurmA=$SLURM_ARRAY_TASK_ID 'NR==slurmA {print $0}' ${output_dir}/pop_list.txt)

## Make directory
mkdir -p $output_dir/$pop

## Use Population (waterbody + ecotype)
grep -w -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv | 
   awk -F ',' -v OFS='\t' -v pop=$pop 'NR!=1 && $10==pop { print $1, $10 }' > $output_dir/$pop/${pop}_pop_file.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/$pop/${pop}_pop_file.txt > $output_dir/$pop/${pop}_ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env
# Filter to those specific samples
bcftools view -S $output_dir/$pop/${pop}_ind_file.txt $vcf | \
    bcftools view --min-ac 1:minor -O z -o $output_dir/$pop/$pop.vcf.gz

# Deactivate bcftools enviroment
conda deactivate

conda activate easySFS-env
python ~/apps/easySFS/easySFS.py -a -i $output_dir/$pop/$pop.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj.txt

# Convert to table
tail -n 2 $output_dir/$pop/${pop}_proj.txt | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop}_proj_long.txt

# Extract best projection number
topproj=$(awk '{print $2}' $output_dir/$pop/${pop}_proj_long.txt | sort -n | tail -1)
bestproj=$(awk -v topproj=$topproj '$2==topproj {print $1 }' $output_dir/$pop/${pop}_proj_long.txt | sort -n | tail -1)

# Run easySFS
python ~/apps/easySFS/easySFS.py -i $output_dir/$pop/$pop.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt -a -f -o $output_dir/$pop/SFS/ --prefix ${pop} --proj $bestproj