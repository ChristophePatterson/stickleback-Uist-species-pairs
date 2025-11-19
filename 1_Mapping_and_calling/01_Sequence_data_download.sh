#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --time=100:00:00
#SBATCH --job-name=sequence_data_download
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# possible sharepoint sites
# MacCollLab1
# MacCollLab2
# OrgOne

# load the Rclone module
module load rclone-uon/1.65.2
# Check rclone remotes
rclone listremotes

# Create output in shared folder
output_dir=(~/data/sticklebacks/seq/)
mkdir -p $output_dir
mkdir -p $output_dir/raw_download
mkdir -p $output_dir/seq_data

cd $output_dir

## Create list of files from 
bigdata=(/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv)
awk -F "," '{print $5 "/" $2}' $bigdata > seq_list.txt
awk -F "," '{print $5 "/" $3}' $bigdata >> seq_list.txt

## Filter to specific samples
# awk -F ',' '$9=="OLAV" {print $0}' $bigdata | awk -F "," '{print $5 "/" $2}' > seq_list.txt
# awk -F ',' '$9=="OLAV" {print $0}' $bigdata | awk -F "," '{print $5 "/" $3}' >> seq_list.txt

# Remove excess file path and split into the two sharepoint back up location (has to be done separately)
grep "sites/MacCollSticklebackLab/Shared Documents/" seq_list.txt |  sed 's/^.\{45\}//' > MacCollSticklebackLab_seq_files.txt
grep "sites/MacColl_stickleback_lab_2/Shared Documents/" seq_list.txt |  sed 's/^.\{49\}//' > MacColl_stickleback_lab_2_seq_files.txt

## Check which files have been downloaded already and cut to basename only
find $output_dir/seq_data/ -wholename '*.gz' | awk -F "/" '{print $NF}' > existing_seq_files.txt
# Create temporary file to hold list of files to download
grep -v -f existing_seq_files.txt MacCollSticklebackLab_seq_files.txt > MacCollSticklebackLab_seq_files_to_download.tmp.txt
grep -v -f existing_seq_files.txt MacColl_stickleback_lab_2_seq_files.txt > MacColl_stickleback_lab_2_seq_files_to_download.tmp.txt

## Number of files to download
echo "Number of files to download from MacCollSticklebackLab:"
wc -l MacCollSticklebackLab_seq_files_to_download.tmp.txt
echo "Number of files to download from MacColl_stickleback_lab_2:"
wc -l MacColl_stickleback_lab_2_seq_files_to_download.tmp.txt

# Create raw download directory
mkdir -p $output_dir/raw_download/

# Copy over wanted list of sequence files
# From MacCollSticklebackLab
rclone --bwlimit 100M --checkers 4 --transfers 4 --onedrive-chunk-size 5M copy --files-from MacCollSticklebackLab_seq_files_to_download.tmp.txt MacCollSticklebackLab: $output_dir
echo "completed download for MacCollSticklebackLab"

# From MacColl_stickleback_lab_2
rclone --bwlimit 100M --checkers 4 --transfers 4 --onedrive-chunk-size 5M copy --files-from MacColl_stickleback_lab_2_seq_files_to_download.tmp.txt MacCollSticklebackLab_2: $output_dir/raw_download/
echo "completed download for MacColl_stickleback_lab_2"

## Find all downloaded files and transfer to single directory
find $output_dir/raw_download/ -name "*.gz" -exec mv {} $output_dir/seq_data/ \;

# unload the rclone module
module unload rclone-uon/1.65.2
