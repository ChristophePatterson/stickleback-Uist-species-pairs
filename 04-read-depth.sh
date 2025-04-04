#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-10
#SBATCH --mem=5g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_readdepth
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0

# extract the individual name variable from sample name files
# Data on all samples
# Define the bigdata file
bigdata="/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-03-28.csv"

# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

# Extract individual using awk
individual=$(awk -F ',' 'BEGIN { OFS="," } { gsub(/^ *| *$/, "", $1); if (FNR == ENVIRON["SLURM_ARRAY_TASK_ID"]) print $1 }' $bigdata)

# Check the result
echo "Individual extracted: $individual"

# set variables
in_filepath=(~/data/sticklebacks/bams/raw_bams)
out_filepath=(~/data/sticklebacks/bams/bamstats/QC/raw_bams)
mkdir -p $out_filepath

echo "This is job $SLURM_ARRAY_TASK_ID and will use sample $individual using bam read $in_filepath/${individual}_raw.bam"

## Run qualimap
~/apps/qualimap_v2.3/qualimap bamqc -bam $in_filepath/${individual}_raw.bam \
    -nt $SLURM_ARRAY_TASK_ID -outformat HTML -outdir $out_filepath/${individual}/ --outfile ${individual}_raw.html


######## # calculate depth for all bams
######## samtools depth \
######## -a \
######## -J \
######## -H \
######## $in_filepath/${individual}_raw.bam |
######## awk -F '\t' '(NR==1) {split($0,header);N=0.0;next;} {N++;for(i=3;i<=NF;i++) a[i]+=int($i);} END { for(x in a) print header[x], a[x]/N;}' > ~/data/sticklebacks/bams/bamstats/${individual}_raw_coverage_depth.txt
######## 
######## # Then in the console run
######## # copy all the depth statistics to a single file
######## ## cat ~/data/sticklebacks/bams/bamstats/*_raw_coverage_depth.txt > ~/data/sticklebacks/bams/bamstats/raw_bam_coverage_depth_all.txt
######## 
######## # and get rid of the file extension and path leaving just the individual name and the mean depth
######## ## sed -i 's/\.bam//' ~/data/sticklebacks/bams/bamstats/raw_bam_coverage_depth_all.txt
######## ## sed -i 's@.*/@@' ~/data/sticklebacks/bams/bamstats/raw_bam_coverage_depth_all.txt
######## ## sed -i 's/_raw//' ~/data/sticklebacks/bams/bamstats/raw_bam_coverage_depth_all.txt
######## 
######## # unload samtools
######## module unload samtools-uoneasy/1.18-GCC-12.3.0