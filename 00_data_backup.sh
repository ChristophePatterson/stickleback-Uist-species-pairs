#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --time=01:00:00
#SBATCH --job-name=StickleBackup
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# possible sharepoint sites
# MacCollLab1
# MacCollLab2
# OrgOne

# load the Rclone module
module load rclone-uon/1.65.2

# Copy all of the files from your folder on Ada to a folder on sharepoint
rclone --transfers 1 --checkers 1 --bwlimit 100M --checksum copy ~/data/sticklebacks MacColl_stickleback_lab_2:Christophe/data/sticklebacks

## Copy over code
rclone --transfers 1 --checkers 1 --bwlimit 100M --checksum copy ~/code MacColl_stickleback_lab_2:Christophe/code


# and check that the two folders are identical
rclone check --one-way ~/data/sticklebacks MacColl_stickleback_lab_2:Christophe/data/sticklebacks

# unload the rclone module
module unload rclone-uon/1.65.2

