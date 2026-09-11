#!/bin/bash
# Laura Dean and Christophe Patterson
# 31/01/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30g
#SBATCH --time=18:00:00
#SBATCH --job-name=vcf2GenomicsGeneral
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
# load conda enviroment that contains bcftools
source /gpfs01/home/${USER}/.bashrc

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
call_pars=(ploidy_aware_HWEPops_MQ10_BQ20)
vcf_ver=($genome_name/$call_pars)

## Sex chromosome ID
Xchr="CM102094.1"
Ychr=""
miti="CM102118.1"

# Activate genomics general
conda activate genomics-general-p3.13

## Convert vcf to geno (all samples and all sites)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz

## Convert vcf to geno (all samples)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz

## Convert vcf to geno (for all samples with mask for coding regions)
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.vcf.gz \
--skipIndels --threads $SLURM_CPUS_PER_TASK | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.masked.geno.gz

## Convert X,Y, and PAR to geno format
# PAR
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz \
--skipIndels --include $Xchr --threads $SLURM_CPUS_PER_TASK --ploidy 2 | bgzip > $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz
# X
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz \
    --skipIndels --include $Xchr --threads $SLURM_CPUS_PER_TASK --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt --ploidyMismatchToMissing | \
    bgzip > $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz

# PAR
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.vcf.gz \
--skipIndels --include $Xchr --threads $SLURM_CPUS_PER_TASK --ploidy 2 | bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz
# X
python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.vcf.gz \
    --skipIndels --include $Xchr --threads $SLURM_CPUS_PER_TASK --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt --ploidyMismatchToMissing | \
    bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz
# Y
## python ~/apps/genomics_general/VCF_processing/parseVCFs.py -i $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.vcf.gz \
##     --skipIndels --include NC_053233.1 --threads $SLURM_CPUS_PER_TASK --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_Y.txt --ploidyMismatchToMissing| \
##     bgzip > $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.Y.geno.gz
## 