#!/bin/bash
# Christophe Patterson
# 20/03/25
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40g
#SBATCH --time=01:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# Module load
module load fastqc-uoneasy/0.12.1-Java-11
module load rclone-uon/1.65.2

#  Input and output directory
dir_input=(~/data/sticklebacks/seq)
dir_output=(~/data/sticklebacks/fastqc)

# Create output directory
mkdir -p $dir_output

# Run fastqc
fastqc -o $dir_output -t 32 $dir_input/*.fastq.gz

# Activate conda
source /gpfs01/home/${USER}/.bashrc
conda activate stickleback-pairs-Uist

# Make specific output folder
mkdir -p $dir_output/multiqc/
# Run multiqc
multiqc -o $dir_output/multiqc/ $dir_output/*fastqc.zip --interactive --pdf

## Push multiqc files to sharepoint (currently html report breaks so need to remove check sum from this rclone)
rclone --ignore-checksum --ignore-size --transfers 1 --checkers 1 --bwlimit 100M copy $dir_output/multiqc MacColl_stickleback_lab_2:Christophe/data/sticklebacks/multiqc/

conda deactivate