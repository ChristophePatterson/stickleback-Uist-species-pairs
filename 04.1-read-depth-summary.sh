#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_readdepth_sum
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.3.3-gfbf-2023b-rstudio 

# set variables
in_filepath=(~/data/sticklebacks/bams/raw_bams)
out_filepath=(~/data/sticklebacks/bams/bamstats/QC/raw_bams)
mkdir -p $out_filepath

# Define the bigdata file
bigdata="/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-03-28.csv"

## Create input config for qualimap
## Reduce total number samples going in using head command (remove if not needed)
awk -F ',' -v filepath="$out_filepath" '{ print $1 " " filepath "/" $1 "/" }' $bigdata | head -n 600 > qualimap.tmp.txt

~/apps/qualimap_v2.3/qualimap multi-bamqc -d qualimap.tmp.txt \
    -outformat HTML -outdir $out_filepath/${individual}/ -outfile global_raw_report.html


