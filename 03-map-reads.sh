#!/bin/bash
# Christophe Patterson & Laura Dean
# 21/03/25
# for running on the UoN HPC Ada

# use devq to test code use defq for running jobs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --array=1-117
#SBATCH --mem=80g
#SBATCH --time=72:00:00
#SBATCH --job-name=bwa_mapping
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load the necessary modules
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0

########################
# Output and input directory

input_directory=(~/data/sticklebacks/seq)
dir_output=(~/data/sticklebacks/bams/raw_bams)
# Make output directory
mkdir -p $dir_output

# Draft genome to use
genome=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna)

# index genome only needs doing once
## bwa index $genome

# Data on all samples
bigdata=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-03-28.csv)
pairdata=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv)

# Create paired data if not already made (cant run for each array as errors arise when files are written at same time)
if [ ! -f $pairdata ]; then
    # All samples from paired reads
	grep -E "DUIN|OBSE|LUIB|CLAC" $bigdata > $pairdata
    # Good coverage outgroups
    grep -E "mara22044|alm222070|NOVSC043|NOVSC116|QbcSGH220044|QbcSGH220048|Lubec003|Lubec001|Ice22054|Ice22041" $bigdata >> $pairdata
fi

# Gets specific sample to work with on this array
individual=$(awk -F ',' 'BEGIN { OFS="," } { gsub(/^ *| *$/, "", $1); if (FNR == ENVIRON["SLURM_ARRAY_TASK_ID"]) print $1 }' $pairdata)
# individual=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $pairdata | awk -F ',' '{ print $1 }')
forward_read=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $pairdata | awk -F ',' '{ print $5 "/" $2 }')
backward_read=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $pairdata | awk -F ',' '{ print $5 "/" $3 }')

# Use awk to process the path
forward_read=$(echo "$forward_read" | awk '{sub(/^sites\/MacColl_stickleback_lab_2\/Shared Documents\//, ""); sub(/^sites\/MacCollSticklebackLab\/Shared Documents\//, ""); print}')
backward_read=$(echo "$backward_read" | awk '{sub(/^sites\/MacColl_stickleback_lab_2\/Shared Documents\//, ""); sub(/^sites\/MacCollSticklebackLab\/Shared Documents\//, ""); print}')

echo "This is job $SLURM_ARRAY_TASK_ID and should use sample $individual"
echo "And should use seq files $input_directory/${forward_read} and $input_directory/${backward_read}" 
echo "Making bam file $dir_output/${individual}_raw.bam"

## Test if file has already been created
## Once created remove raw sequence files
if test -f "$dir_output/${individual}_raw.bam.bai"; then
    echo "${individual} already completed."
    scancel "$SLURM_JOB_ID"
else
    echo "${individual} not mapped: running bwa."
fi

# Extract the header line of the fastq file
file_info=$(zcat $input_directory/${forward_read} | head -n 1)
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
    -v 2 \
    -R "@RG\tID:"$ID"\tSM:"$individual"\tPL:ILLUMINA\tLB:"$individual"\tPU:"$PU"" \
    $genome \
    $input_directory/${forward_read} \
    $input_directory/${backward_read} |
    # mark the duplicate reads then
    # Sort the SAM files
    samtools fixmate --threads $SLURM_CPUS_PER_TASK -m -O BAM - - |
    samtools sort --threads $SLURM_CPUS_PER_TASK -o $dir_output/${individual}_raw.bam

# Index bam file
samtools index $dir_output/${individual}_raw.bam

# Generate info  - look at how well the reads mapped
echo "The raw reads for $individual mapped with the following success:"
samtools flagstat --threads $SLURM_CPUS_PER_TASK $dir_output/${individual}_raw.bam

## Once created remove raw sequence files
if test -f "$dir_output/${individual}_raw.bam.bai"; then
    echo "${individual} successful."
else
    echo echo "${individual} unsuccessful."
fi

cd ~