#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20g
#SBATCH --time=02:00:00
#SBATCH --job-name=VolcanoFinder
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/volcanofinder-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
source ~/.bashrc
module purge

conda activate genomics-general-p3.13

############################
   # SET UP VARIBLES #
############################

### set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback

## Create output directory
output_dir=($wkdir/results/VolcanoFinder)
mkdir -p $output_dir

pop_level=("Population")

## Combine all varibles into single run name
run_name=$(echo "${species}.VolcanoFinder_lv${pop_level}")

############################
   # Create input file #
############################

## Test is Geno file has been created (may need to run whole of 10-sliding window code)
if [ -f $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz ]; then
   echo "Geno.gz file already exists"
else
   echo "Geno.gz file does not exists. Cancelling run."
   scancel "$SLURM_JOB_ID"
fi

### Loop through every population and create allele frequecing
## Create list of individuals to use (# include | shuf | head -n 10) to reduce sample input. Do include iceland for the moment)
if [[ $pop_level == "Population" ]]; then
   grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
      grep -v -E 'Iceland' | awk -F ',' -v OFS='\t' '{ print $1, $10}' > $output_dir/pop_file_${run_name}.txt
else
if [[ $pop_level == "Ecotype" ]]; then
   grep -f $wkdir/vcfs/${species}_subset_samples_withOG.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | \
      grep -v -E 'Iceland' | awk -F ',' -v OFS='\t' '{ print $1, $13}' > $output_dir/pop_file_${run_name}.txt
   fi
fi

# All unique populations
awk '{ print $2 }' $output_dir/pop_file_${run_name}.txt | sort | uniq > $output_dir/uniquePops_${run_name}.txt

## For each unique population create allele frequency file
awk '{ print $2 }' $output_dir/pop_file_${run_name}.txt | sort | uniq | while read Pop
do
   echo "Calculating allele freq for $Pop"
   mkdir -p $output_dir/${Pop}
   ## Run freq.py for each population
   python ~/apps/genomics_general/freq.py -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
   -p $Pop --popsFile $output_dir/pop_file_${run_name}.txt --target minor --asCounts --genoFormat phased --minData 1 \
   --threads $SLURM_CPUS_PER_TASK -o $output_dir/${Pop}/${run_name}_${Pop}.freq

   ## Calculate how many samples are from specic population
   PopN=$(awk '{ print $2 }' $output_dir/pop_file_${run_name}.txt | grep ${Pop}  | wc -l)

   ## Convert pop files in 
   cat $output_dir/${Pop}/${run_name}_${Pop}.freq | \
   awk -v popn="$PopN" 'BEGIN {OFS="\t"} NR==1 {print "scaffold", "position", "x", "n", "folded"} NR!=1 { print $1, $2, $3, popn * 2 , 1}' > $output_dir/${Pop}/${run_name}_${Pop}.VFfreq
   
   ## Create individual VFfreq files for each scafold
   awk 'NR!=1 { print $1 }' $output_dir/${Pop}/${run_name}_${Pop}.VFfreq | sort | uniq | while read scaf
   do
      ## 
      echo "Creating VFFreq for ${Pop} and $scaf"
      grep "${scaf}" $output_dir/${Pop}/${run_name}_${Pop}.VFfreq | awk 'BEGIN {OFS="\t"} NR==1 {print "position", "x", "n", "folded"} NR!=1 { print $2, $3, $4, $5 }' > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq
         
      ## Calucation Empirical unnormalized site frequency spectrum
      # Get total number of sites (need to minus one in awk command to avoid including header)
      nSNPs=$(wc -l  $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq | awk '{ print $1 }')
      cat $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq | awk 'NR!=1 {print $2}' | sort -n | uniq -c | \
      awk -v nsnps=$nSNPs 'BEGIN {OFS="\t"} { print $2, $1/(nsnps-1) }' > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.sfs


   ## Remove unneeded files
   #rm $output_dir/${Pop}/${run_name}_${Pop}.freq
   #rm $output_dir/${Pop}/${run_name}_${Pop}.VFfreq
done








