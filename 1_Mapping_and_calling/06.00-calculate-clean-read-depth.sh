#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-147
#SBATCH --mem=62g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_clean_readdepth
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0
## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.3.3-gfbf-2023b-rstudio 

# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

# extract the individual name variable from sample name files
# Data on all samples
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
# Input path
in_filepath=(~/data/sticklebacks/bams/$genome_name/clean_bams)

# Define the pairdata file
# Input bam files
bam_list="/gpfs01/home/mbzcp2/data/sticklebacks/bams/$genome_name/bamstats/QC/raw_bams/Multi-Bam-QC/HiQ_bam_files.txt"
# extract the individual name variable from sample name files
individual=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $bam_list | awk '{ print $1 }')
bam_file=($in_filepath/$individual.bam)

# Check the result
echo "Individual extracted: $individual"
echo "Bamfile to use is: $bam_file"

# set output variable

out_filepath=(~/data/sticklebacks/bams/$genome_name/bamstats/QC/clean_bams)
mkdir -p $out_filepath

echo "This is job $SLURM_ARRAY_TASK_ID and will use sample $individual using bam read ${bam_file}"

## Run qualimap
~/apps/qualimap_v2.3/qualimap bamqc -bam ${bam_file} \
    -nt $SLURM_ARRAY_TASK_ID --java-mem-size=61G -outdir $out_filepath/${individual}/ -outformat HTML -outfile "${individual}"