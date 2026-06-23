#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25g
#SBATCH --time=24:00:00
#SBATCH --job-name=chrI-inv
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out


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
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

########################
  # ChrI inversion  # 
########################

output_dir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/Regions_of_interest/ChrI_Inv
mkdir -p $output_dir

# Subset bcf to region surronding chrI (add in Uist22521,Uist22542 if you want outgroups)
bcftools view -V indels -r CM102076.1:26500000-27130000 \
  -s Uist22617,Uist22631,Uist22628,Uist22616,Uist22627,Uist22629,Uist22618,Uist22619,Uist22620,Uist22635,Uist22625,Uist22632,Uist22521,Uist22542 \
  $wkdir/vcfs/$vcf_ver/stickleback.bcf | \
  bcftools +fill-tags -- -t AN,AC,AF,MAF |
  bcftools view -Oz -o ${output_dir}/stickleback_DUIN_chrI_inv.vcf.gz
tabix ${output_dir}/stickleback_DUIN_chrI_inv.vcf.gz

bcftools view --min-ac 1:minor -i 'N_ALT>=1' -Oz -o ${output_dir}/stickleback_DUIN_chrI_inv_SNPs.vcf.gz ${output_dir}/stickleback_DUIN_chrI_inv.vcf.gz
tabix ${output_dir}/stickleback_DUIN_chrI_inv_SNPs.vcf.gz

# Highest coverage DUIN Uist22631
# Highest coverage DUIM Uist22617

bcftools view -s Uist22631,Uist22617 -Oz -o ${output_dir}/stickleback_DUIN_chrI_inv_SNPs_top2.vcf.gz ${output_dir}/stickleback_DUIN_chrI_inv_SNPs.vcf.gz
tabix ${output_dir}/stickleback_DUIN_chrI_inv_SNPs_top2.vcf.gz

# Extract atp1a1a coding region
# extract regions contained with atp1a1a (g1117 in DUKE genome and LOC120823921 in Darwin)
grep "g1118.t1" /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf | \
    awk -v OFS='\t' '$3=="CDS" {print $1, $4, $5}' > ${output_dir}/atp1a1a_CDS.bed

# Subset VCF to CDS region
bcftools view -R ${output_dir}/atp1a1a_CDS.bed -Ov -o ${output_dir}/stickleback_DUIN_chrI_inv_SNPs_atp1a1a.vcf.gz ${output_dir}/stickleback_DUIN_chrI_inv_SNPs.vcf.gz 

# Run sliding window analysis
Rscript ~/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/19.1-ChrI-inv-divergence-est.R \
    ${output_dir}/stickleback_DUIN_chrI_inv.vcf.gz \
    ${output_dir}/stickleback_DUIN_chrI_inv_SNPs_atp1a1a.vcf.gz
