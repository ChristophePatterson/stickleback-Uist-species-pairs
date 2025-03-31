#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --time=72:00:00
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

## Create list of files from 
bigdata=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-03-28.csv)
awk -F "," '{print $5 "/" $2}' $bigdata > seq_list.txt
awk -F "," '{print $5 "/" $3}' $bigdata >> seq_list.txt

# Remove excess file path and split into the two sharepoint back up location (has to be done separately)
grep "sites/MacColl_stickleback_lab_2/Shared Documents/" seq_list.txt |  sed 's/^.\{49\}//' > MacColl_stickleback_lab_2_seq_files.txt
grep "sites/MacCollSticklebackLab/Shared Documents/" seq_list.txt |  sed 's/^.\{45\}//' > MacCollSticklebackLab_seq_files.txt

# Copy over wanted list of sequence files
# From MacCollSticklebackLab
rclone --bwlimit 100M --checkers 4 --transfers 4 --onedrive-chunk-size 5M -q copy --files-from MacCollSticklebackLab_seq_files.txt MacCollSticklebackLab: $output_dir
# From MacColl_stickleback_lab_2
rclone --bwlimit 100M --checkers 4 --transfers 4 --onedrive-chunk-size 5M -q copy --files-from MacColl_stickleback_lab_2_seq_files.txt MacColl_stickleback_lab_2: $output_dir

# unload the rclone module
module unload rclone-uon/1.65.2
