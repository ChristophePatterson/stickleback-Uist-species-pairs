#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=01:00:00
#SBATCH --job-name=sliding-window
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/twisst-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
source ~/.bashrc
# load modules
## Remove modules to remove python conflict
module purge
# Careful that the packages you load dont load python
module load R-uoneasy/4.3.2-gfbf-2023a
module load PhyML/3.3.20220408-foss-2023a

## Activate twisst-ete3-p3-6
echo "activate twisst-ete3-p3-6"
conda activate twisst-ete3-p3-6

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback

## Create output directory
output_dir=($wkdir/results/twisst)
mkdir -p $output_dir
########################
## Choose window width
########################
mywindow=(100)

## Test is Geno file has been created (may need to run whole of 10-sliding window code)
if [ -f $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz ]; then
   echo "Geno.gz file already exists"
else
   echo "Geno.gz file does not exists. Cancelling run."
   scancel "$SLURM_JOB_ID"
fi

## Create list of individuals to use (# include | shuf | head -n 10) to reduce sample input. Do include iceland for the moment)
grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
   grep -v -E 'Iceland' | awk -F ',' -v OFS='\t' '{ print $1, $13}' > $output_dir/pop_file_w$mywindow.txt
# Include Iceland
grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
   grep -E 'Iceland' | awk -F ',' -v OFS='\t' '{ print $1, $13}' | sed 's/NA/Ice/' >> $output_dir/pop_file_w$mywindow.txt

awk '{ print $1 }' $output_dir/pop_file_w$mywindow.txt > $output_dir/ind_file_w$mywindow.txt

## Add in new phased sample names 
awk -v OFS='\t' '{ print $1"_A", $2}' $output_dir/pop_file_w$mywindow.txt > $output_dir/phased_pop_file_w$mywindow.txt
awk -v OFS='\t' '{ print $1"_B", $2}' $output_dir/pop_file_w$mywindow.txt >> $output_dir/phased_pop_file_w$mywindow.txt

## Create individual file with just sample name
awk '{ print $1 }' $output_dir/phased_pop_file_w$mywindow.txt > $output_dir/phased_ind_file_w$mywindow.txt

cat $output_dir/pop_file_w$mywindow.txt

### Run Genomics general script for calculating trees over a sliding window
## Best to use one thread because it doesnt take that line
## NOTE THIS CODE HAS BEEN MOVED INTO THE MAIN GENOMICS GENERAL DIRECTORY SO IT CAN ACCESS THE GENOMICS.PY SCRIPT
python ~/apps/genomics_general/phyml_sliding_windows.py -T 1 -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.geno.gz \
  --prefix $output_dir/${species}.phyml_bionj.w$mywindow --windType sites --model GTR --windSize $mywindow -O 0 -M 1 -Ms 1 --indFile $output_dir/ind_file_w$mywindow.txt

## Running Twisst (require install of ete3)
python ~/apps/twisst/twisst.py -t $output_dir/${species}.phyml_bionj.w$mywindow.trees.gz -w $output_dir/${species}.phyml_bionj.w$mywindow.weights.tsv.gz \
   -g anad -g resi -g fw -g Ice --method complete --groupsFile $output_dir/phased_pop_file_w$mywindow.txt


module purge
## Load R module (DO NOT MOVE TO START OF SCRIPT AS IT BREAKS THE PYTHON VERSION)
module load R-uoneasy/4.2.1-foss-2022a
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/11.1-twisst.R $output_dir/${species}.phyml_bionj.w$mywindow