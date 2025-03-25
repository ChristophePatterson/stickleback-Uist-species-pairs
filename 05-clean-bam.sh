#!/bin/bash
# Laura Dean and Christophe
# 22/11/24
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --array=1-165%5
#SBATCH --mem=35g
#SBATCH --time=02:00:00
#SBATCH --job-name=BD_clean
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

########################################################
# SET UP YOUR ENVIRONMENT AND SPECIFY DATA INFORMATION #
########################################################

# 1-10%5 runs 10 arrays with ID 1 to 10 but limits the total number of running at the same time to 5 using the %

# extract the individual name variable from sample name files
individual=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample_names.txt)

# set the input data location
master_filepath=(~/data/sticklebacks/bams)

echo "This is array task ${SLURM_ARRAY_TASK_ID}, cleaning individual $individual, cleaned output BAM files will be written to the folder $master_filepath/clean_bams"

## Test if file has already been created
## Once created remove raw sequence files
if test -f "$master_filepath/clean_bams/${individual}.bam.bai"; then
    echo "${individual} already completed."
    rm -f $master_filepath/raw_bams/${individual}_raw.bam
    rm -f $master_filepath/raw_bams/${individual}_raw.bam.bai
    scancel "$SLURM_JOB_ID"
else
    echo "${individual} not mapped: running bwa."
fi

# load the necessary modules
module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0

# make the output directory if it doesn't already exist
mkdir -p $master_filepath/clean_bams

# Remove unmapped reads and do quality filtering
# -q mapping quality greater than or equal to 40
# -f include reads mapped in a propper pair
# -F Only include reads which are not read unmapped or mate unmapped
samtools view \
--threads $SLURM_CPUS_PER_TASK \
-q 40 \
-f 2 \
-F 4 \
-b $master_filepath/raw_bams/${individual}_raw.bam |
# Mark duplicate reads
samtools markdup -r --threads $SLURM_CPUS_PER_TASK - $master_filepath/clean_bams/$individual.bam
# adding the -r flag to the command above will remove the duplicate reads

# index the final BAM files
samtools index -@ $SLURM_CPUS_PER_TASK $master_filepath/clean_bams/$individual.bam

# Test is work was completed
if test -f "$master_filepath/clean_bams/${individual}.bam.bai"; then
    echo "${individual} completed."
    rm -f $master_filepath/raw_bams/${individual}_raw.bam
    rm -f $master_filepath/raw_bams/${individual}_raw.bam.bai
else
    echo "${individual} failed to clean bam file."
fi

# check the mapping
echo "after cleaning and filtering the final mapping success was:"
samtools flagstat $master_filepath/clean_bams/$individual.bam