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
vcf_ver=ploidy_aware_PopHWE

# load conda enviroment that contains bcftools
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env

#########################
# Concatenate vcf files #
#########################

# write a list of files to be concatenated

if [ ! -f $wkdir/vcfs/$vcf_ver/${species}_ChrLevelVcfFileList.txt ]; then
   ls $wkdir/vcfs/$vcf_ver/${species}_NC*_sorted.bcf > $wkdir/vcfs/$vcf_ver/${species}_ChrLevelVcfFileList.txt
fi
# Concatenate individual chromosome level VCF files
bcftools concat \
    --file-list $wkdir/vcfs/$vcf_ver/${species}_ChrLevelVcfFileList.txt \
    -o $wkdir/vcfs/$vcf_ver/${species}.bcf \
    -O b \
    --threads $SLURM_CPUS_PER_TASK

# index the concatenated VCF file
bcftools index $wkdir/vcfs/$vcf_ver/${species}.bcf
tabix $wkdir/vcfs/$vcf_ver/${species}.bcf

#######################
#### SNP Filtering ####
#######################

bcftools view -v snps -t ^NC_053230.1,NC_053233.1 $wkdir/vcfs/$vcf_ver/${species}.bcf | \
    # Mark GT with less than DP 5 as missing
    bcftools filter -S . -e 'FMT/DP<5' | \
        # Remove SNPs that have average DP of less than 5, greater DP  than 200 and a quality score of less than 60
    bcftools view -e 'AVG(FMT/DP)<5 || AVG(FMT/DP)>200 || QUAL<60' | \
    # Remove SNPs that are missing is more than 80% of samples
    bcftools view -e 'AN/2<N_SAMPLES*0.8' | \
    # Remove SNPs that have a minor allele frequency of less than 2
    bcftools view --min-ac 2:minor -O b -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf

##  Convert to vcf
bcftools view -O z -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

##### Radomly sample one SNP per 1000bp window for rapid assesment of filtering
echo '6. SNPS randomly thinned to one per 10000 bases'
bcftools view -v snps $wkdir/vcfs/$vcf_ver/${species}.bcf | bcftools +prune -n 1 -N rand -w 10000bp -O v -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.rand10000.vcf.gz

##########################
##### LD calculation #####
##########################
# Requires SNP library to have been created

conda deactivate
#  Load specific modules
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig
module load R-uoneasy/4.2.1-foss-2022a

# Go to new directory
cd $wkdir/vcfs/$vcf_ver/

mkdir ld_decay
# move into it
cd ld_decay

# calc ld with plink
plink --vcf $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz --double-id --allow-extra-chr \
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
        $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf -O b -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf

# Convert to vcf
bcftools view -O z $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz

##### Radomly sample one SNP per 10000bp window
echo '6. SNPS randomly thinned to one per 10000 bases'
bcftools +prune -n 1 --nsites-per-win-mode rand -w 10000bp $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf -Ob -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand10000.vcf.gz

##### Radomly sample one SNP per 1000bp window (for masked SNPs)
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.bcf -Oz -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.rand1000.vcf.gz

##### Radomly sample one SNP per 1000bp window
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf -Oz -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz


### Create input for genomics general
## Remove non-species pair samples
# Get list of samples from file
bcftools query -l $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz > $wkdir/vcfs/$vcf_ver/${species}_samples.txt 
# Get info from species pairs data and filter to only include those that are from correct waterbodies
grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
    grep -E 'DUIN|OBSE|LUIB|CLAC' | \
    awk -F ',' '{ print $1 } ' > $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt 

## Filter out non species pair samples and remove sites that are nolonger variable for non masked samples
bcftools view -S $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.bcf | \
    bcftools view --min-ac 2:minor -O z -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz
    # Index
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz

##
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.vcf.gz -Oz -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.rand1000.vcf.gz

## Get geno with outgroup (using Iceland Samples)
grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
    grep -E 'DUIN|OBSE|LUIB|CLAC|OLAV|Iceland|Lubec' | \
    awk -F ',' '{ print $1 } ' > $wkdir/vcfs/$vcf_ver/${species}_subset_samples_withOG.txt 

## Not masked
# Filter to those specific samples
bcftools view -S $wkdir/vcfs/$vcf_ver/${species}_subset_samples_withOG.txt $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz | \
    bcftools view --min-ac 2:minor -O z -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz
# Index
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz


#########################################
 #### Create vcf files for XY & Par ####
#########################################

## Create list of males
sex_der=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/Genomic_sex_determination.txt)
awk '$2=="M" {print $1}' $sex_der | grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt > $wkdir/vcfs/$vcf_ver/male_samples.txt
awk '$2=="F" {print $1}' $sex_der | grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt > $wkdir/vcfs/$vcf_ver/female_samples.txt
## Create ploidy for X (non-PAR region) and Y for input into genomics general
awk -v OFS='\t' '$2=="M" {print $1, 1} $2=="F" {print $1, 2}' $sex_der | grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt > $wkdir/vcfs/$vcf_ver/ploidy_X.txt
awk -v OFS='\t' '$2=="M" {print $1, 1} $2=="F" {print $1, 0}' $sex_der | grep -f $wkdir/vcfs/$vcf_ver/male_samples.txt > $wkdir/vcfs/$vcf_ver/ploidy_Y.txt

## Sum ploidy to get filter params for X chromosome
X_AN=$(awk -F'\t' '{sum+=$2;} END{print sum;}' $wkdir/vcfs/$vcf_ver/ploidy_X.txt)
Y_AN=$(awk -F'\t' '{sum+=$2;} END{print sum;}' $wkdir/vcfs/$vcf_ver/ploidy_Y.txt)

## Create vcf for just PAR, X, and Y (lowers coverage threshold compared to autosomes)
## Further filter to
# PAR
bcftools view -v snps -r NC_053230.1:1-2500000 $wkdir/vcfs/$vcf_ver/${species}.bcf | \
    # Mark GT with less than DP 5 as missing
    bcftools filter -S . -e 'FMT/DP<5' | \
    # Remove SNPs that have average DP of less than 5, greater DP  than 200 and a quality score of less than 60
    bcftools view -e 'AVG(FMT/DP)<5 || AVG(FMT/DP)>200 || QUAL<60' | \
    bcftools view -e 'AN/2<N_SAMPLES*0.8' | \
    bcftools view --min-ac 2:minor -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz
# X
bcftools view -v snps -r NC_053230.1:2500001-20580295 $wkdir/vcfs/$vcf_ver/${species}.bcf | \
    bcftools filter -S . -e 'FMT/DP<2' | \
    bcftools view -e 'AVG(FMT/DP)<2 || AVG(FMT/DP)>200 || QUAL<60' | \
    bcftools view -e "AN<${X_AN}*0.8" | \
    bcftools view --min-ac 2:minor -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz
## Create vcf for just Y
bcftools view -S $wkdir/vcfs/$vcf_ver/male_samples.txt -v snps -r NC_053233.1 $wkdir/vcfs/$vcf_ver/${species}.bcf | \
    bcftools filter -S . -e 'FMT/DP<2 | GT=="het"'| \
    bcftools view -e 'AVG(FMT/DP)<2 || AVG(FMT/DP)>200 || QUAL<60' | \
    bcftools view -e "AN<${Y_AN}*0.8" | \
    bcftools view --min-ac 2:minor -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.vcf.gz
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.vcf.gz

## Merge Autosome and X chromsome calls
# Create merge file list
echo -e "$wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz\n$wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz\n$wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz" \
    > $wkdir/vcfs/$vcf_ver/auto_plus_X_vcf_list.txt
## Merge using bcftools concat
bcftools concat \
    --file-list $wkdir/vcfs/$vcf_ver/auto_plus_X_vcf_list.txt \
    -o $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX.vcf.gz \
    -O z \
    --threads $SLURM_CPUS_PER_TASK
tabix $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX.vcf.gz



#Deactivate env
conda deactivate
### Convert vcfs into genomics general geno formmat
#  Activea genomics general env
conda activate genomics-general-p3.13

## Convert vcf to geno (all samples)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz

## Convert vcf to geno (for all samples with mask for coding regions)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.geno.gz

## Convert X,Y, and PAR to geno format
# PAR
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz \
--skipIndels --include NC_053230.1 --threads $SLURM_CPUS_PER_TASK --ploidy 2 | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz
# X
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz \
    --skipIndels --include NC_053230.1 --threads $SLURM_CPUS_PER_TASK --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt --ploidyMismatchToMissing | \
    bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz
# Y
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.vcf.gz \
    --skipIndels --include NC_053233.1 --threads $SLURM_CPUS_PER_TASK --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_Y.txt --ploidyMismatchToMissing| \
    bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.geno.gz


