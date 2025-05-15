#!/bin/bash
# Christophe Patterson
# 14/05/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50g
#SBATCH --time=0-23:00:00
#SBATCH --job-name=Stickle_call_TTmethod
#SBATCH --array=1-10
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# export the paths and load required environment modules
module load bcftools-uoneasy/1.18-GCC-13.2.0

reference_genome=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna)
# Caculate sequence length of genome (only needs doing once in console)
## bioawk -c fastx '{ print $name, length($seq) }' $reference_genome > $reference_genome.seq_len.txt
# Only keep chromosomes with lengths greater than X
## awk '$2>1000000' $reference_genome.seq_len.txt | awk '{print $1}' > $reference_genome.chrom_names.txt
## wc -l $reference_genome.chrom_names.txt

## Set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback
top_cov=($wkdir/results/TTmethod/vcfs/top_coverage_samples.txt)

## Get name of sample to be called and its corrisponding bamfile
sample=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $top_cov | awk '{ print $1 }')
bamfile=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $top_cov | awk '{ print $2 }')

## Create folder for that specific sample
mkdir -p $wkdir/results/TTmethod/vcfs/$sample

chr_num=0
cat $reference_genome.chrom_names.txt | while read chr 
do
    ## Add one to chromsome number
    chr_num=$((chr_num+1))
    # Print to check everythings working
    echo "Mapping sample $sample to scaffold $chr which is chromosome $chr_num"
    ## Call all sites for this sample
    bcftools mpileup \
    --threads $SLURM_CPUS_PER_TASK \
    --output-type u \
    --max-depth 10000 \
    --min-MQ 20 \
    --min-BQ 30 \
    --platforms ILLUMINA \
    --annotate FORMAT/DP,FORMAT/AD \
    --fasta-ref $reference_genome \
    --regions $chr \
        $bamfile |
    # -m = use the multiallelic caller
    # -v = output variant sites only (TTmethod requires all sites to be called)
    # -P = mutation rate (set at default)
    # -G - = Ignore HWE when making calls
    bcftools call \
    --threads $SLURM_CPUS_PER_TASK \
    -m \
    -P 1e-6 \
    -a GQ \
    -O b \
    -G - \
    -o $wkdir/results/TTmethod/vcfs/${sample}/${sample}_chr${chr_num}.bcf

    # Sort and Index the BCF file
    bcftools sort -O b -o $wkdir/results/TTmethod/vcfs/${sample}/${sample}_sorted_chr${chr_num}.bcf $wkdir/results/TTmethod/vcfs/${sample}/${sample}_chr${chr_num}.bcf
    bcftools index $wkdir/results/TTmethod/vcfs/${sample}/${sample}_sorted_chr${chr_num}.bcf

    ## Remove unsorted bcf if sorted file have been created
    if [ ! -f $wkdir/results/TTmethod/vcfs/${sample}/${sample}_sorted_chr${chr_num}.bcf.csi ]; then
    	rm $wkdir/results/TTmethod/vcfs/${sample}/${sample}_chr${chr_num}.bcf
    fi

    ## Convert to vcf
    bcftools view -O z $wkdir/results/TTmethod/vcfs/${sample}/${sample}_sorted_chr${chr_num}.bcf > $wkdir/results/TTmethod/vcfs/${sample}/${sample}_chr${chr_num}.vcf.gz
    tabix $wkdir/results/TTmethod/vcfs/${sample}/${sample}_chr${chr_num}.vcf.gz
    
done

