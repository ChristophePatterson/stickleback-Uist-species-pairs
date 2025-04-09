#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-117
#SBATCH --mem=62g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_readdepth
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0
## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.3.3-gfbf-2023b-rstudio 

# extract the individual name variable from sample name files
# Data on all samples
# Define the pairdata file
pairdata="/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv"

# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

# Extract individual using awk
individual=$(awk -F ',' 'BEGIN { OFS="," } { gsub(/^ *| *$/, "", $1); if (FNR == ENVIRON["SLURM_ARRAY_TASK_ID"]) print $1 }' $pairdata)

# Check the result
echo "Individual extracted: $individual"

# set variables
in_filepath=(~/data/sticklebacks/bams/raw_bams)
out_filepath=(~/data/sticklebacks/bams/bamstats/QC/raw_bams)
mkdir -p $out_filepath

echo "This is job $SLURM_ARRAY_TASK_ID and will use sample $individual using bam read $in_filepath/${individual}_raw.bam"

## Run qualimap
~/apps/qualimap_v2.3/qualimap bamqc -bam $in_filepath/${individual}_raw.bam \
    -nt $SLURM_ARRAY_TASK_ID --java-mem-size=61G -outdir $out_filepath/${individual}/ -outformat HTML -outfile "${individual}_raw"