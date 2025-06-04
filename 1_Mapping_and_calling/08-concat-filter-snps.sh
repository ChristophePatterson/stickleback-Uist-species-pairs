#!/bin/bash
# Laura Dean and Christophe Patterson
# 31/01/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=20g
#SBATCH --time=18:00:00
#SBATCH --job-name=stickle_cat_vcfs
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# load conda enviroment that contains bcftools
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env

#########################
# Concatenate vcf files #
#########################

# write a list of files to be concatenated

if [ ! -f $wkdir/${species}_ChrLevelVcfFileList.txt ]; then
   ls $wkdir/vcfs/${species}_NC*_sorted.bcf > $wkdir/${species}_ChrLevelVcfFileList.txt
fi
# Concatenate individual chromosome level VCF files
bcftools concat \
--file-list $wkdir/${species}_ChrLevelVcfFileList.txt \
-o $wkdir/vcfs/${species}.bcf \
-O b \
--threads $SLURM_CPUS_PER_TASK

# index the concatenated VCF file
bcftools index $wkdir/vcfs/${species}.bcf
tabix $wkdir/vcfs/${species}.bcf

#####################
# VARIANT FILTERING #
#####################

# Extract only the SNPs from the VCF file (as it contains indels as well)
bcftools view -v snps $wkdir/vcfs/${species}.bcf -O b -o $wkdir/vcfs/${species}_SNPs.bcf

# Index the SNP only VCF file
bcftools index $wkdir/vcfs/${species}_SNPs.bcf
tabix $wkdir/vcfs/${species}_SNPs.bcf


#######################
#### SNP Filtering ####
#######################

##### Radomly sample one SNP per 1000bp window for rapid assesment of filtering
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 10000bp $wkdir/vcfs/${species}_SNPs.bcf -Ov -o $wkdir/vcfs/${species}_SNPs.rand10000.vcf.gz

# All Varients
echo '0. All varients'
bcftools view $wkdir/vcfs/${species}_SNPs.bcf | grep -v -c '^#'

echo "Number of samples"
bcftools query -l $wkdir/vcfs/${species}_SNPs.bcf | wc -l

# All snps
echo '1. All single nucleotides'
bcftools view $wkdir/vcfs/${species}_SNPs.bcf | grep -v -c '^#'

# SNPs genotyped in more than 

# Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60

echo '2. Depth greater than 5, average depth for each sample greater than 10 and less than 200, quality score greater than 60'
bcftools filter -S . -e 'FMT/DP<5' $wkdir/vcfs/${species}_SNPs.bcf | \
bcftools view -e 'AVG(FMT/DP)<5 || AVG(FMT/DP)>200 || QUAL<60' -O b > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.bcf | grep -v -c '^#'

##### Radomly sample one SNP per 1000bp window for rapid assesment of filtering
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 10000bp $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.rand10000.bcf
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.rand10000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.rand10000.vcf.gz

# Removes SNPs that are not present in more than 80% samples
echo "4. Removing SNPs that arn't genotyped in more than 80% samples"
bcftools view -e 'AN/2<N_SAMPLES*0.8' -O b $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.bcf | grep -v -c '^#'

# Removing SNPs with a minor allele freq less than 0.05
echo "5. Removing alleles that have a minor allele count of less that 2"
bcftools view --min-ac 2 -O b $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz


##########################
##### LD calculation #####
##########################
# Requires SNP library to have been created

conda deactivate
#  Load specific modules
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig
module load R-uoneasy/4.2.1-foss-2022a

# Go to new directory
cd $wkdir/vcfs/

mkdir ld_decay
# move into it
cd ld_decay

# calc ld with plink
plink --vcf $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
-r2 gz --ld-window 1000 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out ${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.LD

# Calculate the cororlation over distance 
# use chmod u+x to make sure file can be excisuted 
/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/LD_dist_calc_python3.py -i ${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.LD.ld.gz -o ${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.LD

# Plot results using R
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/LD_dist_calc_plot.R ${species} ${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.LD.ld_decay_bins
# Clear modules
module purge

################################################
    ### Mask (remove coding regions) ###
##################################################

# Reactivate conda env
conda activate bcftools-env

# Create masking file for input into -T agugment
##awk ' NR!=1 {OFS="\t";print $1, $2, $3}' /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GAculeatus_UGA_version5_genomic_annotations.tsv > /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GAculeatus_UGA_version5_genomic_annotations.bcftools_mask.txt
## Remove all SNPs that are contained with the regions specificied in the file
bcftools view -T ^/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GAculeatus_UGA_version5_genomic_annotations.bcftools_mask.txt \
        $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf -O b -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf

# Convert to vcf
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz
tabix -f $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz

# Record how many records are on each file
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf | grep -v -c '^#'
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf | grep -v -c '^#'

cd $wkdir/vcfs/

##### Radomly sample one SNP per 10000bp window
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 --nsites-per-win-mode rand -w 10000bp $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand10000.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand10000.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand10000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand10000.vcf.gz

##### Radomly sample one SNP per 1000bp window
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.bcf
bcftools view -H $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.vcf.gz

##### Radomly sample one SNP per 100000bp window
echo '6. SNPS randomly thinned to one per 1000000 bases'
bcftools +prune -n 1 -N rand -w 100000bp $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand100000.bcf
bcftools view -H $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand100000.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand100000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand100000.vcf.gz

### Create input for genomics general
## Remove non-species pair samples
# Get list of samples from file
bcftools query -l $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz > $wkdir/vcfs/${species}_samples.txt 
# Get info from species pairs data and filter to only include those that are from correct waterbodies
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
    grep -E 'DUIN|OBSE|LUIB|CLAC' | \
    awk -F ',' '{ print $1 } ' > $wkdir/vcfs/${species}_subset_samples.txt 

## Filter out non species pair samples and remove sites that are nolonger variable
bcftools view -S $wkdir/vcfs/${species}_subset_samples.txt $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz | \
    bcftools view --min-ac 2[minor] -O z -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.vcf.gz
# Index
tabix $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.vcf.gz

## Random SNPs for just paired populations
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.vcf.gz -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.rand1000.bcf
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.rand1000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.rand1000.vcf.gz

## Filter out non species pair samples and remove sites that are nolonger variable
bcftools view -S $wkdir/vcfs/${species}_subset_samples.txt $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf | \
    bcftools view --min-ac 2[minor] -O z -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz
    # Index
tabix $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.9.MAF2_SpPair.vcf.gz

## Get geno with outgroup (using Iceland Samples)
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
    grep -E 'DUIN|OBSE|LUIB|CLAC|OLAV|Iceland|Lubec' | \
    awk -F ',' '{ print $1 } ' > $wkdir/vcfs/${species}_subset_samples_withOG.txt 

# Filter to those specific samples
bcftools view -S $wkdir/vcfs/${species}_subset_samples_withOG.txt $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz | \
    bcftools view --min-ac 2[minor] -O z -o $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair-wOG.vcf.gz
# Index
tabix $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair-wOG.vcf.gz

#Deactivate env
conda deactivate
### Convert vcfs into genomics general geno formmat
#  Activea genomics general env
conda activate genomics-general-p3.13

## Convert vcf to geno (for just species pairs - masked)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair.geno.gz

## Convert vcf to geno (for just species pairs - retaining coded regions)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz

## Convert vcf to geno (for just species pairs but with outgroup)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair-wOG.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked_SpPair-wOG.geno.gz

## Convert vcf to geno (for all samples)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.geno.gz
