#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=12:00:00
#SBATCH --array=1-3
#SBATCH --job-name=sliding-window-pca
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## Due to overlap in writing files, if slurm array is not equal to 1 then wait 15 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware_HWEPops_MQ10_BQ20

########################
  # LEA - PCA & SNMF  # 
########################

## Run variables
vcf=${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX
vcf_full=$wkdir/vcfs/$vcf_ver/$vcf.vcf.gz
wndsize=25000
wndslid=5000

output_dir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/$vcf
mkdir -p $output_dir


if [ $SLURM_ARRAY_TASK_ID == "1" ]; then
   ## Get list of chromosomes to use
   bcftools query -f '%CHROM\n' $vcf_full | sort | uniq > $output_dir/chrom_list.txt 
   ### Populations to use
   bcftools query -l $vcf_full > $output_dir/samples_in_vcf.txt
   echo -e "OBSE\nDUIN\nLUIB\nCLAC" > $output_dir/Pops_interest.txt
   #### Extract sample information
   awk -F ',' -v OFS='\t' '{ print $1, $13 ,$9}' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
         grep -f $output_dir/samples_in_vcf.txt | \
         grep -f $output_dir/Pops_interest.txt > $output_dir/pop_file.txt
   awk '{print $1}' $output_dir/pop_file.txt > $output_dir/samples.txt
fi

# Get chromosome number for array
chr=$(awk "FNR==$SLURM_ARRAY_TASK_ID" $output_dir/chrom_list.txt)
## Create output director
mkdir -p $output_dir/$chr

# Subset to specific chromosome
bcftools view -r $chr -S $output_dir/samples.txt --min-ac 2:minor -O z -o $output_dir/$chr/stickleback.$chr.vcf.gz $vcf_full

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
        $output_dir/$chr/stickleback.$chr.vcf.gz $vcf_ver $wndsize $wndslid "TRUE"

# Remove subset vcf
rm $output_dir/$chr/stickleback.$chr.vcf.gz
        
## Analysis where population can be changed (specifically removing CLAC)
## run_name=(noCLAC)
## output_dir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/$vcf.$run_name
## mkdir -p $output_dir
## 
## # Get names of samples in vcf
## bcftools query -l $wkdir/vcfs/$vcf_ver/$vcf.vcf.gz > $output_dir/samples_in_vcf.txt
## # Create file of waterbodies that are of interes
## echo -e "OBSE\nDUIN\nLUIB" > $output_dir/Pops_interest.txt
## # Get sample information and subset to those water bodies
## awk -F ',' -v OFS='\t' '{ print $1, $13 ,$9}' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
##          grep -f $output_dir/samples_in_vcf.txt | \
##          grep -f $output_dir/Pops_interest.txt > $output_dir/pop_file.txt
## # Get subset list of sample files
## awk '{print $1}' $output_dir/pop_file.txt > $output_dir/samples.txt
## 
## # Remove CLAC samples and remove nolonger variable snps
## bcftools view -S $output_dir/samples.txt --min-ac 2:minor -O z -o $output_dir/$vcf.$run_name.vcf.gz $wkdir/vcfs/$vcf_ver/$vcf.vcf.gz
## 
## # Run analysis again
## Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
##         $output_dir/$vcf.$run_name.vcf.gz $vcf_ver $wndsize $wndslid "TRUE"
## 
## rm $output_dir/$vcf.$run_name.vcf.gz


