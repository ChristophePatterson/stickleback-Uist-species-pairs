#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=18:00:00
#SBATCH --job-name=sliding-window
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module load PhyML/3.3.20220408-foss-2023a
module load R-uoneasy/4.2.1-foss-2022a

# Open conda enviroment
conda activate twisst-stickleback

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

## Create output directory
mkdir -p $wkdir/results/twisst

## Test is Geno file has been created (may need to run whole of 10-sliding window code)
if [ -f $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz ]; then
   echo "Geno.gz file already exists"
else
   echo "Geno.gz file does not exists. Cancelling run."
   scancel "$SLURM_JOB_ID"
fi

## Create list of individuals to use (# include | shuf | head -n 10 to reduce sample input and not include iceland for the moment)
grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
   grep -v -E 'Iceland' | shuf | head -n 20 |  awk -F ',' -v OFS='\t' '{ print $1, $13}' > $wkdir/results/twisst/pop_file.txt
# Include Iceland
grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
   grep -E 'Iceland' | awk -F ',' -v OFS='\t' '{ print $1, $13}' | sed 's/NA/Ice/' >> $wkdir/results/twisst/pop_file.txt

## Add in new phased sample names 
awk -v OFS='\t' '{ print $1"_A", $2}' $wkdir/results/twisst/pop_file.txt > $wkdir/results/twisst/phased_pop_file.txt
awk -v OFS='\t' '{ print $1"_B", $2}' $wkdir/results/twisst/pop_file.txt >> $wkdir/results/twisst/phased_pop_file.txt

## Create individual file with just sample name
awk '{ print $1 }' $wkdir/results/twisst/pop_file.txt > $wkdir/results/twisst/ind_file.txt

cat $wkdir/results/twisst/pop_file.txt

### Run Genomics general script for calculating trees over a sliding window
## Best to use one thread because it doesnt take that line
## NOTE THIS CODE HAS BEEN MOVED INTO THE MAIN GENOMICS GENERAL DIRECTORY SO IT CAN ACCESS THE GENOMICS.PY SCRIPT
python ~/apps/genomics_general/phyml_sliding_windows.py -T 1 -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.geno.gz \
   --prefix $wkdir/results/twisst/${species}.phyml_bionj.w50 --windType sites --model GTR --windSize 1000 -O 0 -M 1 -Ms 1 --indFile $wkdir/results/twisst/ind_file.txt

### Running Twisst (require install of ete3)
python ~/apps/twisst/twisst.py --help

python ~/apps/twisst/twisst.py -t $wkdir/results/twisst/${species}.phyml_bionj.w50.trees.gz -w $wkdir/results/twisst/${species}.weights.csv.gz -o $wkdir/results/twisst/${species}.topologies.trees\
   -g anad -g resi -g fw -g Ice --method complete --groupsFile $wkdir/results/twisst/phased_pop_file.txt
