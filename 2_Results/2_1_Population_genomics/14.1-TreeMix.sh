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
#SBATCH --job-name=TreeMix
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out


## Load modules 
module load boost-uoneasy/1.82.0-GCC-13.2.0
module load gsl-uoneasy/2.7-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load vcftools-uoneasy/0.1.16-GCC-12.3.0
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig

## Run test tree mix

inputdir=(/gpfs01/home/mbzcp2/data/sticklebacks/vcfs) # input file locations
outputdir=(/gpfs01/home/mbzcp2/data/sticklebacks/results/treemix) # output file locations
# Create output location
mkdir -p $outputdir

## Conversion of vcf to treemix input file

# staring vcf
vcf=(stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair)

## Step 1 remove all missing data
bcftools view -e 'AN/2<N_SAMPLES' -O z $inputdir/$vcf.vcf.gz > $inputdir/$vcf.noN.vcf.gz
bcftools view $inputdir/$vcf.noN.vcf.gz | grep -v -c '^#'

## Step 2 filter to remove LD
bcftools +prune -n 1 -N rand -w 1000bp $inputdir/$vcf.noN.vcf.gz > $outputdir/$vcf.noN.r1000.vcf.gz
bcftools view $inputdir/$vcf.noN.r1000.vcf.gz | grep -v -c '^#'

## Step 3 create clust file
grep -f $inputdir/stickleback_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
    awk -F ',' '{ print $1, $1, $10}' > $outputdir/pop_file.clust

## Convert from vcf to tree mix file using code provide from 
# https://speciationgenomics.github.io/Treemix/ # Need to direct where plink to tree mix can be found
# & https://bitbucket.org/nygcresearch/treemix/downloads/


conda activate treemix-p2
/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/vcf2treemix.sh $outputdir/$vcf.noN.r1000.vcf.gz $outputdir/pop_file.clust
conda deactivate

### Run treemix

module load R-uoneasy/4.2.1-foss-2022a
for i in {1..5} 
do
    treemix -i $outputdir/$vcf.noN.r1000.treemix.frq.gz -m $i -o $outputdir/stickleback.$i
    Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/14.2-TreeMix-plot.R $outputdir/stickleback.$i
done