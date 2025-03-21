#!/bin/bash
# Christophe Patterson & Laura Dean
# 21/03/25
# for running on the UoN HPC Ada

# use devq to test code use defq for running jobs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=19
#SBATCH --array=6-152
#SBATCH --mem=35g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_mapping
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load the necessary modules
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0

########################
# Output and input directory

input_directory=(/share/MacColl_shared/Christophe/seq)
dir_output=(/share/MacColl_shared/Christophe/bam/raw_bams)
# Make output directory
mkdir -p $dir_output

# Draft genome to use
genome=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna)

# index genome only needs doing once
## bwa index $genome

# Gets specific sample to work with on this array
individual=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample_names.txt)
echo "This is job $SLURM_ARRAY_TASK_ID and should use sample $individual"
echo "And should use seq files $input_directory/${individual}_R1.fastq.gz and $input_directory/${individual}_R2.fastq.gz" 

# Extract the header line of the fastq file
file_info=$(zcat $input_directory/${individual}_R1.fastq.gz | head -n 1)
# Save the pieces of information you need as variables
flowcell_ID=$(cut -d ":" -f3 <<< "$file_info")
lane_no=$(cut -d ":" -f4 <<< "$file_info")
sample_barcode=$(cut -d ":" -f10 <<< "$file_info")

# store the read group information for the ID and PU fields as variables from the individual ones you just created
PU=$flowcell_ID.$lane_no.$sample_barcode
ID=$flowcell_ID.$lane_no

bwa mem \
    -t $SLURM_CPUS_PER_TASK \
    -M \
    -R "@RG\tID:"$ID"\tSM:"$individual"\tPL:ILLUMINA\tLB:"$individual"\tPU:"$PU"" \
    $genome \
    $input_directory/${individual}_R1.fastq.gz \
    $input_directory/${individual}_R2.fastq.gz |
    # mark the duplicate reads then
    # Sort the SAM files
    samtools fixmate --threads $SLURM_CPUS_PER_TASK -m -O BAM - - |
    samtools sort --threads $SLURM_CPUS_PER_TASK -o $dir_output/${individual}_raw.bam

# Index bam file
samtools index $dir_output/${individual}_raw.bam

# Generate info  - look at how well the reads mapped
echo "The raw reads for $individual mapped with the following success:"
samtools flagstat --threads 19 $dir_output/${individual}_raw.bam

cd ~