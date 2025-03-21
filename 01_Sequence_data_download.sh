#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --time=48:00:00
#SBATCH --job-name=sequence_data_download
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# possible sharepoint sites
# MacCollLab1
# MacCollLab2
# OrgOne

# load the Rclone module
module load rclone-uon/1.65.2

# Create output in shared folder
output_dir=(~/data/sticklebacks/seq/)
mkdir -p $output_dir

# Get complete list of all files in back up folder
rclone lsf MacColl_stickleback_lab_2:HPC_data_backup/bigdata/trimmed_fqs > seq_list.txt
## Just get sequence files with Uist in name
grep Uist seq_list.txt > seq_list_uist.txt
awk -F "\"*,\"*" '$8==2' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_metadata_ADCM_MB.csv | \awk -F "\"*,\"*" {'print $2'} > sample_names.txt
awk '{print $1 "_R1.fastq.gz"}' sample_names.txt > sample_pairs_seq.txt
awk '{print $1 "_R2.fastq.gz"}' sample_names.txt >> sample_pairs_seq.txt

# Copy over wanted list of sequence files
rclone --bwlimit 100M --checkers 4 --transfers 4 --onedrive-chunk-size 5M -q copy --files-from sample_pairs_seq.txt MacColl_stickleback_lab_2:HPC_data_backup/bigdata/trimmed_fqs $output_dir

# unload the rclone module
module unload rclone-uon/1.65.2
