#!/bin/bash
# Laura Dean and Christophe
# 22/11/24
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=19
#SBATCH --array=1-5
#SBATCH --mem=35g
#SBATCH --time=02:00:00
#SBATCH --job-name=BD_clean
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

########################################################
# SET UP YOUR ENVIRONMENT AND SPECIFY DATA INFORMATION #
########################################################

# extract the individual name variable from sample name files
individual=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample_names.txt)

# set the input data location
master_filepath=(~/data/sticklebacks/bams)

echo "This is array task ${SLURM_ARRAY_TASK_ID}, cleaning individual $individual, cleaned output BAM files will be written to the folder $master_filepath/cleaned_bams"

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

# check the mapping
echo "after cleaning and filtering the final mapping success was:"
samtools flagstat ~/data/sticklebacks/bams/$individual.bam