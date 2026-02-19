#!/bin/bash
# Christophe Patterson
# 18/02/26
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20g
#SBATCH --time=23:00:00
#SBATCH --job-name=Stickle_annotate
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

########################################################
# SET UP YOUR ENVIRONMENT AND SPECIFY DATA INFORMATION #
########################################################


## Extract CDS regions from the GTF file
## load conda enviroment that contains gffread
source /gpfs01/home/${USER}/.bashrc
conda activate gffread-env

genome_dir=/gpfs01/home/mbzcp2/data/sticklebacks/genomes

# Set variables
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
gtf_file=$genome_dir/Duke_GAcu_1.gtf
reference_genome=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/$genome_name)

# Output directory
output_dir=$genome_dir/${genome_name}_functional_annotation
mkdir -p $output_dir # create output directory

## Because reference genome and gtf have different chromosome names, we need to change the chromosome names in the gtf file to match those in the reference genome
# Get sequencing report for the reference genome
seq_report=$genome_dir/${genome_name}_sequence_report.tsv

# Replace the chromosome names in the gtf file with those in the reference genome
# First filter out "chr" from chromosome names in the gtf file, a nd remove leading "0" from chromosome names in the gtf file
awk '{ gsub("chr", "", $1); gsub(/^0+/, "", $1); print }' $gtf_file > $output_dir/${genome_name}_ChrNames.gtf
awk '{print $7}' $seq_report # Print chromosome names from the sequence report to check they match those in the gtf file
awk '{print $4}' $seq_report # Print chromosome names from the sequence report to check they match those in the gtf file
awk '{print $1}' $output_dir/${genome_name}_ChrNames.gtf | uniq # Print chromosome names from the gtf file to check they match those in the sequence report

# For each line in the gtf file, match column 1 of the gtf file with column 4 of the sequence report, and replace it with column 7 of the sequence report
awk -v OFS=',' 'NR==FNR{a[$4]=$7; next} {if($1 in a) $1=a[$1]}1' $seq_report $output_dir/${genome_name}_ChrNames.gtf  > $output_dir/${genome_name}_ChrNames_fixed.gtf
# Replace the commas in the gtf file with tabs, and remove the "gene_id" from the attributes column, and replace it with "gene_id " (with a space after it) to make it compatible with gffread
awk '{ gsub("_id,", "_id "); print }' $output_dir/${genome_name}_ChrNames_fixed.gtf |
    awk '{ gsub(";,gene_id", "; gene_id"); print }' |
    awk '{ gsub(",", "\t"); print }' > $output_dir/${genome_name}_ChrNames_fixed_notab.gtf

# Check that the chromosome names in the gtf file now match those in the sequence report
awk '{print $1}' $output_dir/${genome_name}_ChrNames_fixed_notab.gtf | uniq
grep ">" $reference_genome.fna | awk '{print $1}' | sed 's/>//g' # Check that the chromosome names in the reference genome match those in the gtf file

# head file for testing
# head -n 100 $output_dir/${genome_name}_ChrNames_fixed_notab.gtf > $output_dir/${genome_name}_ChrNames_fixed_notab_head.gtf

## Run gffread to extract the CDS regions from the gtf file, and output them in fasta format
gffread -w $output_dir/${genome_name}_CDS.fa -g $reference_genome.fna $output_dir/${genome_name}_ChrNames_fixed_notab.gtf > $output_dir/${genome_name}_gffread.log 2>&1


# Check that the CDS fasta file has been created and contains the expected number of sequences
grep ">" $output_dir/${genome_name}_CDS.fa | wc -l 

## Extract CDS.fa file from two other stickleback genomes
# From fGasAcu3.hap1.1
awk '$3 != "gene"' $genome_dir/GCF_964276395.1/genomic.gtf | gffread -w $genome_dir/GCF_964276395.1/GCF_964276395.1_fGasAcu3.hap1.1_genomic_CDS.fa -g $genome_dir/GCF_964276395.1/GCF_964276395.1_fGasAcu3.hap1.1_genomic.fna > $genome_dir/GCF_964276395.1/gffread_fGasAcu3.log 2>&1


# From v5
gffread -w $genome_dir/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic_CDS.fa -g $genome_dir/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna $genome_dir/GCF_016920845.1/genomic.gff > $genome_dir/GCF_016920845.1/stickleback_v5_ensembl_genes_gffread.log 2>&1

## Unload conda environment
conda deactivate

## Load blast
module load blast-uoneasy/2.14.1-gompi-2023a

# Blast the CDS sequences from the Duke genome against the CDS sequences from the other two stickleback genomes, to find orthologous genes
# Create blast databases for the other two stickleback genomes
makeblastdb -in $genome_dir/GCF_964276395.1/GCF_964276395.1_fGasAcu3.hap1.1_genomic_CDS.fa -dbtype nucl -out $genome_dir/GCF_964276395.1/GCF_964276395.1_fGasAcu3.hap1.1_genomic_CDS
makeblastdb -in $genome_dir/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic_CDS.fa -dbtype nucl -out $genome_dir/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic_CDS 

# Blast the CDS sequences from the Duke genome against the blast databases for the other two stickleback genomes
blastn -query $output_dir/${genome_name}_CDS.fa -db $genome_dir/GCF_964276395.1/GCF_964276395.1_fGasAcu3.hap1.1_genomic_CDS \
    -out $output_dir/${genome_name}_vs_fGasAcu3_blastn.out -evalue 1e-5 -outfmt "6 qaccver saccver pident evalue qcovs ssciname stitle"
blastn -query $output_dir/${genome_name}_CDS.fa -db $genome_dir/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic_CDS \
    -out $output_dir/${genome_name}_vs_v5_blastn.out -evalue 1e-5 -outfmt "6 qaccver saccver pident evalue qcovs ssciname stitle"

# Load R for parsing blast output and finding best hits
module load R-uoneasy/4.2.1-foss-2022a

# Parse blast output and find best hits for each gene in the Duke genome, and add this information to the gtf file
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/DUKE_genome_functionAnnotation_blast_parsing.R



