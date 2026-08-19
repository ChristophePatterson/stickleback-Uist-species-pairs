#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
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

# Subset bcf to region surronding chrI with all samples
bcftools view -V indels -r CM102076.1:26500000-27130000 \
  $wkdir/vcfs/$vcf_ver/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz | \
  bcftools +fill-tags -- -t AN,AC,AF,MAF |
  bcftools view -Oz -o ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz
tabix ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz


# Extract atp1a1a coding region
# extract regions contained with atp1a1a (g1117 in DUKE genome and LOC120823921 in Darwin)
grep -w "g1118" /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf | \
    awk -v OFS='\t' '$5=="CDS" {print $1, $2, $3}' > ${output_dir}/atp1a1a_CDS.bed

grep -w "g1118" /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf | \
    awk -v OFS='\t' '$5=="CDS" && $7=="t1" {print $1, $2, $3}' > ${output_dir}/atp1a1a_CDS_t1.bed

grep -w "g1118" /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf | \
    awk -v OFS='\t' '$5=="CDS" {print $1, $2, $3, $5, $6, $7, $9}' > ${output_dir}/atp1a1a_CDS_info.bed

grep -w "g1118" /gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCA_046562415.1_Duke_GAcu_1.0_genomic_functional_annotation/GCA_046562415.1_Duke_GAcu_1.0_genomic_blast_matches.gtf | \
    awk -v OFS='\t' '{print $1, $2, $3, $5, $6, $7, $9}' > ${output_dir}/atp1a1a_info.bed


# Subset VCF to CDS region
bcftools view -R ${output_dir}/atp1a1a_CDS.bed -Ov -o ${output_dir}/stickleback_DUIN_chrI_inv_SNPs_atp1a1a.vcf.gz ${output_dir}/stickleback_DUIN_chrI_inv_SNPs.vcf.gz 
bcftools view -r CM102076.1:26836909-26867066 -Ov -o ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.vcf.gz ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz

# Run sliding window analysis
Rscript ~/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/19.1-ChrI-inv-divergence-est.R \
    ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_chrI_inv.vcf.gz \
    ${output_dir}/stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a.vcf.gz

conda deactivate

module load raxml-ng-uoneasy/1.2.0-GCC-12.3.0

## Run RAxML
SNP_library="stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_atp1a1a"
phy_file="$SNP_library.phy"

# Change into output directory
cd $output_dir

# Run RAxML
pwd
echo $phy_file

# raxml-ng --check --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --threads $SLURM_CPUS_PER_TASK --prefix GTGTR4_G_ASC_LEWIS
# raxml-ng --parse --msa $wkdir/vcfs/$vcf_ver/$phy_file --model GTGTR4+G+ASC_LEWIS --prefix GTGTR4_G_ASC_LEWIS
raxml-ng --all --msa $phy_file --model GTGTR4+G+ASC_LEWIS --bs-trees 500 --tree pars{10},rand{10} --threads $SLURM_CPUS_PER_TASK --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10

## Map bootstrap values on the ML tree
raxml-ng --support --tree ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bestTree --bs-trees ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bootstraps --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml --threads $SLURM_CPUS_PER_TASK

## Test convergence of bootstrap trees
raxml-ng --bsconverge --bs-trees ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bootstraps --prefix ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml.bsconverge --threads $SLURM_CPUS_PER_TASK --bs-cutoff 0.01


## module purge
## Load R
module load R-uoneasy/4.2.1-foss-2022a

# Run R plotting script
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/16.1-RAxML_plot.R \
    ${SNP_library}_raxml_GTGTR4_G_ASC_LEWIS_BS500_P10R10.raxml


