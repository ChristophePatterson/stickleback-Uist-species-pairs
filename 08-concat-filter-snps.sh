#!/bin/bash
# Laura Dean and Christophe Patterson
# 31/01/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=100:00:00
#SBATCH --job-name=stickle_cat_vcfs
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# load modules
module load bcftools-uoneasy/1.18-GCC-13.2.0

#########################
# Concatenate vcf files #
#########################

# write a list of files to be concatenated (have to escape the ls command for Augusta as I have an alias for it)
ls $wkdir/vcfs/${species}_NC*.vcf.gz > $wkdir/${species}_ChrLevelVcfFileList.txt

# Concatenate individual chromosome level VCF files
#### bcftools concat \
#### --file-list $wkdir/${species}_ChrLevelVcfFileList.txt \
#### -o $wkdir/vcfs/${species}.vcf.gz \
#### -O z \
#### --threads $SLURM_CPUS_PER_TASK
#### 
#### # index the concatenated VCF file
#### bcftools index $wkdir/vcfs/${species}.vcf.gz
#### tabix -p vcf $wkdir/vcfs/${species}.vcf.gz

#####################
# VARIANT FILTERING #
#####################

# Extract only the SNPs from the VCF file (as it contains indels as well)
###### bcftools view -v snps $wkdir/vcfs/${species}.vcf.gz -O z -o $wkdir/vcfs/${species}_SNPs.vcf.gz

# Index the SNP only VCF file
bcftools index $wkdir/vcfs/${species}_SNPs.vcf.gz
tabix -p vcf $wkdir/vcfs/${species}_SNPs.vcf.gz


#######################
#### SNP Filtering ####
#######################

# All Varients
echo '0. All varients'
bcftools view $wkdir/vcfs/${species}.vcf.gz | grep -v -c '^#'

echo "Number of samples"
bcftools query -l $wkdir/vcfs/${species}_SNPs.vcf.gz | wc -l

# All snps
echo '1. All single nucleotides'
bcftools view $wkdir/vcfs/${species}_SNPs.vcf.gz | grep -v -c '^#'

# SNPs genotyped in more than 

# Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60

echo '2. Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60'
bcftools filter -S . -e 'FMT/DP<10' $wkdir/vcfs/${species}_SNPs.vcf.gz | \
bcftools view -e 'AVG(FMT/DP)<10 || AVG(FMT/DP)>200 || QUAL<60' -O b > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.vcf.gz

##### Radomly sample one SNP per 1000bp window for rapid assesment of filtering
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 10000bp $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.rand10000.bcf
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.rand10000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.rand10000.vcf.gz

# Removes SNPs that are not present in more than 80% samples
echo "4. Removing SNPs that arn't genotyped in more than 80% samples"
bcftools view -e 'AN/2<N_SAMPLES*0.8' -O b $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf | grep -v -c '^#'

# Removing SNPs with a minor allele freq less than 0.05
echo "5. Removing alleles that have a minor allele count of less that 2"
bcftools view --min-ac 2 -O b $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf | grep -v -c '^#'

bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.vcf.gz


##########################
##### LD calculation #####
##########################
######  # Requires SNP library to have been created
######  
######  # Go to new directory
######  cd $wkdir/vcfs/
######  
######  mkdir ld_decay
######  # move into it
######  cd ld_decay
######  
######  # bcftools view -O z $output_dir/$BCF_FILE.snps.bcf > $output_dir/$BCF_FILE.snps.vcf.gz
######  
######  # calc ld with plink
######  plink --vcf $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.vcf.gz --double-id --allow-extra-chr \
######  --set-missing-var-ids @:# \
######  --maf 0.01 --geno 0.1 --mind 0.5 \
######  -r2 gz --ld-window 1000 --ld-window-kb 1000 \
######  --ld-window-r2 0 \
######  --make-bed --out $BCF_FILE.LD
######  
######  # Calculate the cororlation over distance 
######  /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/6_LD_dist_calc_python3.py -i $BCF_FILE.LD.ld.gz -o $BCF_FILE.LD
######  
######  # Plot results using R
######  Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/6_LD_dist_calc_plot.R
######  
######  cd $output_dir


##### Radomly sample one SNP per 1000bp window
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 10000bp $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand10000.bcf
bcftools view $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand10000.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand10000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand10000.vcf.gz

##### Radomly sample one SNP per 1000bp window
echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools +prune -n 1 -N rand -w 1000bp $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf -Ob -o $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.bcf
bcftools view -H $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.bcf | grep -v -c '^#'
bcftools view -O z $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.bcf > $wkdir/vcfs/${species}_SNPs.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz

