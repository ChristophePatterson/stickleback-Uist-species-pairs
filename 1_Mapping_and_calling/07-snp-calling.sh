#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20g
#SBATCH --time=23:00:00
#SBATCH --job-name=Stickle_call
#SBATCH --array=1-22
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

########################################################
# SET UP YOUR ENVIRONMENT AND SPECIFY DATA INFORMATION #
########################################################

# load conda enviroment that contains bcftools
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env

reference_genome=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna)
ploidy_file=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/stickleback_ploidy.txt)

# Caculate sequence length of genome (only needs doing once in console)
## bioawk -c fastx '{ print $name, length($seq) }' $reference_genome > $reference_genome.seq_len.txt
# Only keep chromosomes with lengths greater than X
## awk '$2>1000000' $reference_genome.seq_len.txt | awk '{print $1}' > $reference_genome.chrom_names.txt
## wc -l $reference_genome.chrom_names.txt

# specify your array config file that lists the chromosome numbers
# extract the chromosome name from the array config file
chr=$(awk "NR==$SLURM_ARRAY_TASK_ID" $reference_genome.chrom_names.txt)

# Thresholds for Mapping and base quality
MQthres=10
BQthres=20
# set variables
master_filepath=(~/data/sticklebacks) # set the master data location
master_output=($master_filepath/vcfs/ploidy_aware_HWEPops_MQ${MQthres}_BQ${BQthres}) # Set output lociation
mkdir -p $master_output # create output location

VCF=stickleback_${chr} # set the name of the output vcf file
## regionsdir=/gpfs01/home/mbzlld/code_and_scripts/Regions_files/G_aculeatus


# print to the file the array that is being worked on...
echo "This is array task $SLURM_ARRAY_TASK_ID, calling SNPs for chromosome ${chr}, and writing them to the file ${VCF}.bcf"

############################
# SNP and genotype calling #
############################

# create a list of all of the BAM files that we will call into the same variant file (but only if it doesn't exist)
if [ ! -f $master_output/bamlist.txt ]; then
	cat /gpfs01/home/mbzcp2/data/sticklebacks/bams/bamstats/QC/clean_bams/Multi-Bam-QC/HiQ_bam_files.txt | \
	grep -v "Obsm_641" | awk '{print $2 }' | grep -v "Obsm_641" > $master_output/bamlist.txt
fi

## Create Population file
if [ ! -f $master_output/PopFile.txt ]; then
	cat /gpfs01/home/mbzcp2/data/sticklebacks/bams/bamstats/QC/clean_bams/Multi-Bam-QC/HiQ_bam_files.csv | \
	grep -v "Obsm_641" | awk -F ',' -v OFS='\t' 'NR != 1 {print $1, $14}' > $master_output/PopFile.txt
fi

## Genomic sex 
if [ ! -f $master_output/Gsex.ped ]; then
	cat /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/Genomic_sex_determination.ped  | \
	grep -v "Obsm_641" > $master_output/Gsex.ped
fi

# Generate genotype likelihoods for the BAM files using mpileup
# then pipe this to call to generate a BCF file of genetic variants

### ! Removed below code ###
# --regions-file $regionsdir/Chromosome_${SLURM_ARRAY_TASK_ID}.txt
# -m = use the multiallelic caller
# -v = output variant sites only
# -P = mutation rate (set at default)
# -G - = Ignore HWE when making calls
# -G $master_output/PopFile.txt \

bcftools mpileup \
--threads $SLURM_CPUS_PER_TASK \
--output-type u \
--max-depth 10000 \
--min-MQ $MQthres \
--min-BQ $BQthres \
--platforms ILLUMINA \
--annotate FORMAT/DP,FORMAT/AD \
--fasta-ref $reference_genome \
--regions $chr \
--bam-list $master_output/bamlist.txt |
bcftools call \
--threads $SLURM_CPUS_PER_TASK \
-m \
-v \
-P 1e-6 \
-a GQ \
-O b \
-G $master_output/PopFile.txt \
--ploidy M,F \
--ploidy-file $ploidy_file \
--samples-file $master_output/Gsex.ped \
-o $master_output/$VCF.bcf

#unzip the vcf file
#gzip -d $master_output/$VCF

echo "The SNPs have now been called, proced to sorting and indexing"

# Sort and Index the VCF file
bcftools sort -O b -o $master_output/${VCF}_sorted.bcf $master_output/$VCF.bcf
bcftools index $master_output/${VCF}_sorted.bcf

## Remove unsorted bcf if sorted file have been created
if [ -f $master_output/${VCF}_sorted.bcf.csi ]; then
	rm $master_output/$VCF.bcf
fi

# Output some check information on the VCF file you have generated:
# list the sameples contained in the VCF file
echo "These are the individuals in the VCF file:"
bcftools query -l $master_output/${VCF}_sorted.bcf 
# Count all variants in the file
echo "This is the number of variants in the file:"
bcftools view -H $master_output/${VCF}_sorted.bcf | wc -l
# Count all SNPs in the file
echo "This is the number of SNPs in the file:"
bcftools view -H -v snps $master_output/${VCF}_sorted.bcf | wc -l
# Count all indels in the file
echo "This is the number of indels in the file:"
bcftools view -H -v indels $master_output/${VCF}_sorted.bcf | wc -l
