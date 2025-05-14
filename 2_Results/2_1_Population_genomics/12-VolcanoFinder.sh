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

pop_level=("Ecotype")

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

   zcat $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz | awk 'NR!=1 { print $1 }'  | sort | uniq | while read scaf
   do

      ## Filter geno to specific pop and scaf
      python ~/apps/genomics_general/filterGenotypes.py --threads 1 -i $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
      -o $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.geno.gz -p ${Pop} --popsFile $output_dir/pop_file_${run_name}.txt --minCalls 1 --include ${scaf}

      ## Run freq.py for each population
      python ~/apps/genomics_general/freq.py -g $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.geno.gz \
      -p $Pop --popsFile $output_dir/pop_file_${run_name}.txt --target minor --asCounts --genoFormat phased \
      --threads $SLURM_CPUS_PER_TASK -o $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.freq

      ## Calculate how many samples are in the population
      PopN=$(awk '{ print $2 }' $output_dir/pop_file_${run_name}.txt | grep ${Pop}  | wc -l)

      ## Filter posistion
      awk '{ print $2 }' $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.freq > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.pos

      python ~/apps/genomics_general/freq.py -g $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.geno.gz > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.missNo
      
      ## Get geno file, filter out position that are not included in .freq file, then count number of "N/N" in each line
      #zcat $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.geno.gz | cut -f2- | sed -u '1d' | grep -f $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.pos | \
      #grep "N/N" -o -n | awk -F ':' '{ print $1 }'| uniq -c | awk -v popn=$PopN '{ print (popn-$1)*2, 1}' > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.missNo

      ## Convert pop files in 
      cat $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.freq | \
      awk -v popn="$PopN" 'BEGIN {OFS="\t"} NR!=1 { print $2, $3 }' > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq.tmp
      ## File header
      echo -e "position\tx\t\tn\tfolded" > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq
      ## Combine freq with number of genotype counts
      paste $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq.tmp $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.missNo >> $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq

      rm -f $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq.tmp
      ## Create individual VFfreq files for each scafold

      
      ### ## Calucation Empirical unnormalized site frequency spectrum
      ### # Get total number of sites (need to minus one in awk command to avoid including header)
      ### nSNPs=$(wc -l  $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq | awk '{ print $1 }')
      ### cat $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.VFfreq | awk 'NR!=1 {print $2}' | sort -n | uniq -c | \
      ### awk -v nsnps=$nSNPs 'BEGIN {OFS="\t"} { print $2, $1/(nsnps-1) }' > $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.sfs
### 
      ### ## Get min and max x
      ### xmax=$(tail -1 /gpfs01/home/mbzcp2/data/sticklebacks/results/VolcanoFinder/resi/stickleback.VolcanoFinder_lvEcotype_resi_NC_053233.1.sfs | awk '{ print $1 }')
      ### xmin=$(head -1 /gpfs01/home/mbzcp2/data/sticklebacks/results/VolcanoFinder/resi/stickleback.VolcanoFinder_lvEcotype_resi_NC_053233.1.sfs | awk '{ print $1 }')
### 
      ### ## Generating vlookup table
      ### ## EXAMPLE ./VolcanoFinder –p SpectFile D P MODEL nmin nmax xmin xmax LookupPrefix
      ### ~/apps/volcanofinder_v1.0/VolcanoFinder –p $output_dir/${Pop}/${run_name}_${Pop}_${scaf}.sfs D P MODEL $PopN  xmin xmax LookupPrefix
   done

   ## Remove unneeded files
   #rm $output_dir/${Pop}/${run_name}_${Pop}.freq
   #rm $output_dir/${Pop}/${run_name}_${Pop}.VFfreq
done








