#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20g
#SBATCH --array=1-45
#SBATCH --time=18:00:00
#SBATCH --job-name=sliding-window
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module purge 
module load R-uoneasy/4.2.1-foss-2022a

# Activate conda enviroment
conda activate genomics-general-p3.13

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/sliding-window

# Make folder
mkdir -p $wkdir/results/sliding-window/All_Pop_comparison

## Create popfile and unique combination of population (if doesn't exist already)
if [ ! -f $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt ]; then
    grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
        awk -F ',' '{ print $1, $10}' | sed s/OLST/OLAV/ > $wkdir/results/sliding-window/All_Pop_comparison/pop_file.txt

    # Get all unique populations
    awk '{print $2}' $wkdir/results/sliding-window/All_Pop_comparison/pop_file.txt | sort | uniq > $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt
    ## Get all unique combination of populations
    Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/Get_combinations.R $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt
    wc -l $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt
fi
## Get unique combination of waterbodies
pop1=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $1 }' $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt)
pop2=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $2 }' $wkdir/results/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt)

# Create output directory
outputdir=($wkdir/results/sliding-window/All_Pop_comparison/${pop1}_${pop2})
mkdir -p $outputdir

### Run sliding window
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/All_Pop_comparison/pop_file.txt -p ${pop1} -p ${pop2}

## Plot results
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}"
