#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --time=02:00:00
#SBATCH --job-name=sequence_data_download
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# possible sharepoint sites
# MacCollLab1
# MacCollLab2
# OrgOne

# load the Rclone module
module load rclone-uon/1.65.2

# Get complete list of all files in back up folder
rclone lsf MacColl_stickleback_lab_2:HPC_data_backup/bigdata/trimmed_fqs > seq_list.txt
## Just get sequence files with Uist in name
grep Uist seq_list.txt | head -n 20 > seq_list_uist.txt

# Copy over wanted list of sequence files
rclone copy --files-from seq_list_uist.txt MacColl_stickleback_lab_2:HPC_data_backup/bigdata/trimmed_fqs ~/data/sticklebacks/

# unload the rclone module
module unload rclone-uon/1.65.2
