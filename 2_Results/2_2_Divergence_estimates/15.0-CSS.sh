#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80g
#SBATCH --array=1-21
#SBATCH --time=48:00:00
#SBATCH --job-name=CSS
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out


## Due to overlap in writing files, if slurm array is not equal to 1 then wait 15 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

############################
   # PREPARE ENVIRONMENT #
############################
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
conda activate bcftools-env
# Load R
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware_HWEPops_MQ10_BQ20

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

# Parameters for CSS calculation
wndsize=25000
sliding=5000
wdnmthd="basepair" #Unit of window- and stepsizes as number of SNPs (locus) or base pairs (basepair)"
mnSNP=1
mthd=pca
MAF=0.05

## Create input pop file
output_prefix=stickleback.wnd$wndsize.sld$sliding.mnSNP$mnSNP.mth$wdnmthd-$mthd.MAF$MAF.NoCLAC
output_dir=$wkdir/results/$vcf_ver/sliding-window/CSS/$output_prefix
# Create directory
mkdir -p $output_dir

echo "stickleback.wnd$wndsize sld$sliding mnSNP$mnSNP mth$wdnmthd-$mthd MAF$MAF"

## Script output to location of vcf and needs cleaner file name
# Define which vcf to use
vcf=/gpfs01/home/mbzcp2/data/sticklebacks/vcfs/$vcf_ver/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX.vcf.gz

if [ $SLURM_ARRAY_TASK_ID == "1" ]; then
   ## Get list of chromosomes to use
   bcftools query -f '%CHROM\n' $vcf | sort | uniq > $output_dir/chrom_list.txt 
   ### Populations to use
   bcftools query -l $vcf > $output_dir/samples_in_vcf.txt
   echo -e "OBSE\nDUIN\nLUIB" > $output_dir/Pops_interest.txt
   #### Extract sample information
   grep 
   awk -F ',' -v OFS='\t' '{ print $1, $13 ,$9}' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
         grep -f $output_dir/samples_in_vcf.txt | \
         grep -f $output_dir/Pops_interest.txt > $output_dir/pop_file.txt
   awk '{print $1}' $output_dir/pop_file.txt > $output_dir/samples.txt
fi

## Get chromosome
chr=$(awk "FNR==$SLURM_ARRAY_TASK_ID" $output_dir/chrom_list.txt)
## Copy over to outpute folder
## Add -r NC_053212.1:25500000-27000000 to reduce size in test
## Or -R /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/chr_test.tmp.txt for multiple contigs
bcftools view -r $chr -S $output_dir/samples.txt --min-ac 2:minor -O z -o $output_dir/stickleback.$chr.vcf.gz $vcf

## Remove prior results if they exist
if [ -f $output_dir/stickleback.$chr.gds ]; then
   echo "gds file already exists so removing"
   rm $output_dir/stickleback.$chr.gds 
else
   echo "gds file does not exists, running analysis."
fi

# Change to output directy
cd $output_dir

## Loop through and calculate CSS for each chromosome

## Run Rscript for each chromosome
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/CSSm.R \
            $output_dir/stickleback.$chr.vcf.gz pop_file.txt $wndsize $sliding $mnSNP $wdnmthd $mthd $MAF $SLURM_ARRAY_TASK_ID > $output_dir/CSS_log_$chr.txt

## PCACSSm_permutation.R file.vcf file.CSSm.dmat.gz file.CSSm.txt file.grouplist npermutations"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/CSSm_permutations.R \
   stickleback.$chr.vcf.gz \
   stickleback.$chr.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.dmat.gz \
   stickleback.$chr.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.txt pop_file.txt 10000 > CSS_perm_log_$chr.txt

# Remove tempory vcf for specific chromosome
rm $output_dir/stickleback.$chr.vcf.gz
rm $output_dir/stickleback.$chr.gds

## Merge all Perm files together - if there are 21 files already created
permfilesNo=$(ls $output_dir/stickleback.*.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.10000perm.txt | wc -l)
if [ $permfilesNo == 21 ]; then
   echo "All $permfilesNo, perm files created so merging output from all"
   echo -e "chr\tstart\tend\tnsnps\tcss\tpval" > $output_dir/${output_prefix}.CSSm.10000perm.txt
   awk FNR!=1 $output_dir/stickleback.*.${wndsize}${wdnmthd}${sliding}step.window.${mthd}.pop_file.CSSm.10000perm.txt >> $output_dir/${output_prefix}.CSSm.10000perm.txt
   ## Plot in R
   Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/15.1-CSS_plot.R \
      $output_dir ${output_prefix}.CSSm.10000perm.txt
else
   echo "There are only $permfilesNo permutation files so not merging"
fi


