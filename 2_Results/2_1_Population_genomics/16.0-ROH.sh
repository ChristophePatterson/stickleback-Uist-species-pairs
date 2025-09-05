#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH --time=18:00:00
#SBATCH --job-name=ROH
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
species=stickleback
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/ROH
mkdir -p $outdir

# Set SNP library
SNP_lib=${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair
vcf=$wkdir/vcfs/$vcf_ver/$SNP_lib.vcf.gz

## Get sample names 
bcftools query -l $vcf > $outdir/samples_in_vcf.txt

# Get lengths of entire genome
bioawk -c fastx '{ sum += length($seq) } END { print sum }' /gpfs01/home/mbzcp2/data/sticklebacks/genomes/$genome_name.fna > $outdir/$genome_name.length.txt

## Get population samples
grep -w -f $outdir/samples_in_vcf.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
        awk -F ',' -v OFS='\t' '$13!="st" { print $1, $10, $13 } $13=="st" { print $1, $10, "fw" }' > $outdir/pop_file.txt
# Get unique waterbodies
awk '{ print $2 }' $outdir/pop_file.txt | sort | uniq > $outdir/pop_uniq.txt

echo -e "pop\tsample\thetn\tROHsum" > "$outdir/het_counts.txt"
## Calculate ROH for each population
while read line; do
    # Start loop for each population
    pop=$line
    echo Running ROH for $pop
    ## Get samples names
    grep $pop $outdir/pop_file.txt | awk '{ print $1 }' > $outdir/${pop}_samples.txt
    # Filter vcf to just those samples and remove invarient site
    bcftools view -S $outdir/${pop}_samples.txt $vcf | \
        bcftools view --min-ac 1:minor -O z -o $outdir/${pop}.vcf.gz
    # Calculate ROH
    bcftools roh -G30 $outdir/${pop}.vcf.gz -o $outdir/${pop}.out.txt
    # ROH viz
    roh-viz -i $outdir/${pop}.out.txt -v $outdir/${pop}.vcf.gz -o $outdir/${pop}_rmme.html
    # Split data
    grep "RG" $outdir/${pop}.out.txt > $outdir/${pop}.out.RG.txt
    
    # Get genotypes and limit to het
    bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' $outdir/${pop}.vcf.gz | awk '$4=="0/1"' > $outdir/${pop}.genotypes.txt

    # Calculate number and sum of ROH, plus non-ROH het for each sample
    Rscript Helper_scripts/IBrisk_calc.R $outdir/${pop}.genotypes.txt $outdir/${pop}.out.RG.txt 

done < $outdir/pop_uniq.txt

## Combine into a single data set
## ROHs
awk 'NR != 1 { print }' $outdir/*.out.RG.txt | grep -v '#' > $outdir/All.RG.txt
## Individual Stats
cat $outdir/*.out.RG.calcs.txt > $outdir/All.RG.calcs.txt






