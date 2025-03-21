#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=4g
#SBATCH --time=5:00:00
#SBATCH --array=1-5
#SBATCH --job-name=BD_readdepth
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0

# extract the individual name variable from sample name files
individual=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample_names.txt)

# set variables
in_filepath=/share/MacColl_shared/Christophe/bam/raw_bams

# calculate depth for all bams
samtools depth \
-a \
-J \
-H \
$in_filepath/${individual}_raw.bam |
awk -F '\t' '(NR==1) {split($0,header);N=0.0;next;} {N++;for(i=3;i<=NF;i++) a[i]+=int($i);} END { for(x in a) print header[x], a[x]/N;}' > ~/data/sticklebacks/bams/${individual}_raw_coverage_depth.txt

# Then in the console run
# copy all the depth statistics to a single file
## cat ~/data/sticklebacks/bams/*_raw_coverage_depth.txt > ~/data/sticklebacks/bams/raw_bam_coverage_depth_all.txt

# and get rid of the file extension and path leaving just the individual name and the mean depth
## sed -i 's/\.bam//' ~/data/sticklebacks/bams/raw_bam_coverage_depth_all.txt
## sed -i 's@.*/@@' ~/data/sticklebacks/bams/raw_bam_coverage_depth_all.txt
## sed -i 's/_raw//' ~/data/sticklebacks/bams/raw_bam_coverage_depth_all.txt

# unload samtools
module unload samtools-uoneasy/1.18-GCC-12.3.0