#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --array=1-16
#SBATCH --time=18:00:00
#SBATCH --job-name=private-alleles
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module purge 
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load R-uoneasy/4.2.1-foss-2022a

# Activate conda enviroment
conda activate genomics-general-p3.13

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/sliding-window

# Make folder
mkdir -p $wkdir/results/sliding-window/private-alleles

if [ ! -f $wkdir/results/sliding-window/private-alleles/pop_file.txt ]; then
    # Get list of samples from each Uist Population
    grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
        grep -E 'DUIN|OBSE|LUIB|CLAC|OLAV|TORM' | \
        awk -F ',' '{ print $1, $10}' | sed s/OLST/OLAV/ | sed s/TOST/TORM/ > $wkdir/results/sliding-window/private-alleles/pop_file.txt
## Get list of samples from Uist
    awk '{print $1}' $wkdir/results/sliding-window/private-alleles/pop_file.txt > $wkdir/results/sliding-window/private-alleles/sample_file.txt
## Get all other samples but cluster by region not population
    grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
            grep -v -f $wkdir/results/sliding-window/private-alleles/sample_file.txt | \
            awk -F ',' '{ print $1, $7}' >> $wkdir/results/sliding-window/private-alleles/pop_file.txt
## Remake sample file with all samples
    awk '{print $1}' $wkdir/results/sliding-window/private-alleles/pop_file.txt > $wkdir/results/sliding-window/private-alleles/sample_file.txt
fi

awk '{print $2}' $wkdir/results/sliding-window/private-alleles/pop_file.txt | sort | uniq | wc -l

# Get Population to work with from array
pop1=$(awk '{print $2}' $wkdir/results/sliding-window/private-alleles/pop_file.txt | sort | uniq | awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $1 }')

# Create output directory
outputdir=($wkdir/results/sliding-window/private-alleles/${pop1})
mkdir -p $outputdir

## Subset to all population that arn't in pop1
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
        awk -F ',' '{ print $1, $10}' | grep -v "$pop1" | awk '{print $1}' > $outputdir/non_${pop1}_samples.txt
##  All pop1 samples
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
        awk -F ',' '{ print $1, $10}' | grep "$pop1" | awk '{print $1}' > $outputdir/${pop1}_samples.txt

## Get input vcf
vcf=($wkdir/vcfs/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz)

# subset to samples from pop1 that are still variable
bcftools view -S $outputdir/${pop1}_samples.txt $vcf | \
    bcftools view --min-ac 1[minor] -O b -o $outputdir/$pop1.bcf
bcftools index $outputdir/$pop1.bcf

# Subset to all non pop1 samples that are still variable
bcftools view -S $outputdir/non_${pop1}_samples.txt $vcf | \
    bcftools view --min-ac 1[minor] -O b -o $outputdir/n$pop1.bcf
bcftools index $outputdir/n$pop1.bcf

## Calculate stats
bcftools stats -s - $outputdir/$pop1.bcf > $outputdir/$pop1.stats.txt
# Extract just het stats
grep "PSC" $outputdir/$pop1.stats.txt > $outputdir/$pop1.PSC.stats.txt

# How many variable snps are in pop1 that are not seen in any other sampled individual
bcftools isec -n~10 $outputdir/$pop1.bcf $outputdir/n$pop1.bcf > $outputdir/${pop1}_private_alleles.txt
## Get bcftool of private alleles
bcftools view -R $outputdir/${pop1}_private_alleles.txt $outputdir/$pop1.bcf > $outputdir/$pop1.priv.bcf
## Calculate stats
bcftools stats -s - $outputdir/$pop1.priv.bcf > $outputdir/$pop1.priv.stats.txt
# Extract just het stats
grep "PSC" $outputdir/$pop1.priv.stats.txt > $outputdir/$pop1.priv.PSC.stats.txt

## Plot heterozgousity
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.5-Private-alleles-plot.R $outputdir/$pop1

## Number of allels
num_pri=$(cat $outputdir/${pop1}_private_alleles.txt| wc -l )
total_snps=$(bcftools view --min-ac 1[minor] $vcf | wc -l)

## If data base of private allele numbers doesn't exit add to file
touch $wkdir/results/sliding-window/private-alleles/private_allele_number.txt
## Then add info that file
echo -e $pop1 "\t" $num_pri >> $wkdir/results/sliding-window/private-alleles/private_allele_number.txt




