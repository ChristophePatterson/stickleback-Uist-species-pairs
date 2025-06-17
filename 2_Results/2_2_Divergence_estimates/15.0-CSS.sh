#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80g
#SBATCH --time=48:00:00
#SBATCH --job-name=CSS
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
conda activate bcftools-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

# Parameters for CSS calculation
wndsize=25000
sliding=5000
wdnmthd="locus" #Unit of window- and stepsizes as number of SNPs (locus) or base pairs (basepair)"
mnSNP=1
mthd=pca
MAF=0.05

## Create input pop file
output_dir=$wkdir/results/sliding-window/CSS/stickleback.wnd$wndsize.sld$sliding.mnSNP$mnSNP.mth$wdnmthd-$mthd.MAF$MAF
# Create 
mkdir -p $output_dir

echo "stickleback.wnd$wndsize sld$sliding mnSNP$mnSNP mth$wdnmthd-$mthd MAF$MAF"

## Script output to location of vcf and needs cleaner file name
# Define which vcf to use
vcf=/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/ploidy_aware/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz
## Copy over to outpute folder
## Add -r NC_053212.1:25500000-27000000 to reduce size in test
## Or -R /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/chr_test.tmp.txt for multiple contigs
bcftools view -O z -o $output_dir/stickleback.vcf.gz $vcf

# Query names of samples in VCF
bcftools query -l $vcf > $output_dir/samples.txt

## Create pop file for samples
grep -f $output_dir/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $13, $9}' > $output_dir/pop_file.txt

## Deactivate conda
conda deactivate
# Load R
module load R-uoneasy/4.2.1-foss-2022a
## Plot results

## Remove prior results if they exist
if [ -f $output_dir/stickleback.gds ]; then
   echo "gds file already exists so removing"
   rm $output_dir/stickleback.gds 
else
   echo "gds file does not exists, running analysis."
fi

# Change to output directy
cd $output_dir
### CSSm.R file.vcf file.grouplist windowsize stepsize minsnpperwindow [locus|basepair] [pca|mds] minorallelefrequency
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/CSSm.R \
            stickleback.vcf.gz pop_file.txt $wndsize $sliding $mnSNP $wdnmthd $mthd $MAF > $output_dir/CSS_log.txt

## Merge all .CSSm.dmat.gz together
if [ -f $output_dir/stickleback.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.dmat.gz ]; then
   echo "CSSm file already exists so removing"
   rm $output_dir/stickleback.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.dmat.gz
fi
cat $output_dir/*.CSSm.dmat.gz > $output_dir/stickleback.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.dmat.gz


## PCACSSm_permutation.R file.vcf file.CSSm.dmat.gz file.CSSm.txt file.grouplist npermutations"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/CSSm_permutations.R stickleback.vcf.gz \
   stickleback.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.dmat.gz \
   stickleback.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.txt pop_file.txt 10000 > CSS_perm_log.txt

