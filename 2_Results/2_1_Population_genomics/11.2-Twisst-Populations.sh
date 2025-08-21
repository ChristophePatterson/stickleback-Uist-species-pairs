#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6
#SBATCH --mem=20g
#SBATCH --time=48:00:00
#SBATCH --job-name=twisst-waterbody-pairs
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## Due to overlap in writing files, if slurm array is not equal to 1 then wait 30 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

############################
   # PREPARE ENVIRONMENT #
############################
source ~/.bashrc
# load modules
## Remove modules to remove python conflict
module purge
# Careful that the packages you load dont load python
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

## Create output directory
output_dir=($wkdir/results/$vcf_ver/twisst/Population_comparison)
mkdir -p $output_dir
########################
## Choose window width
########################
mywindow=(100)

## Get unique combination of waterbodies
if [ ! -f $output_dir/pop_uniq.txt_combn.txt ]; then
    ## Get each unique Population (but replacing st with fw)
    grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
        awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $13 } $13=="st" { print $1, $9, "fw" }' > $output_dir/pop_file.txt
    # Get unique waterbodies
    awk '{ print $2 }' $output_dir/pop_file.txt | sort | uniq > $output_dir/pop_uniq.txt
    ## Get all combination of waterbodies
    Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/Get_combinations.R $output_dir/pop_uniq.txt 2
    wc -l $output_dir/pop_uniq.txt_combn.txt
fi

pop1=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $1 }' $output_dir/pop_uniq.txt_combn.txt)
pop2=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $2 }' $output_dir/pop_uniq.txt_combn.txt)

## Make output directory
mkdir -p $output_dir/${pop1}_${pop2}

## Create list of samples to be used
grep -E "${pop1}|${pop2}" $output_dir/pop_file.txt | awk '{print $1, $2 "_" $3}' > $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt

# Set weighting method "complete" or "fixed"
weight_method=("complete")

# Create unique run name

## Test is Geno file has been created (may need to run whole of 10-sliding window code)
if [ -f $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz ]; then
   echo "Geno.gz file already exists"
else
   echo "Geno.gz file does not exists. Cancelling run."
   scancel "$SLURM_JOB_ID"
fi

# ID file for all samples
awk '{ print $1 }' $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt > $output_dir/${pop1}_${pop2}/ind_file_${pop1}_${pop2}.txt
# ID file for females -w matches whole word otherwise CLAM1 match CLAM1 and CLAM10
grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt | awk '{ print $1 }' > $output_dir/${pop1}_${pop2}/ind_file_${pop1}_${pop2}.X.txt

## Add in new phased sample names 
awk -v OFS='\t' '{ print $1"_A", $2}' $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt > $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.txt
awk -v OFS='\t' '{ print $1"_B", $2}' $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt  >> $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.txt

## Add in new phased female sample names 
grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt | awk -v OFS='\t' '{ print $1"_A", $2}' > $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.X.txt
grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt $output_dir/${pop1}_${pop2}/${pop1}_${pop2}_popfile.txt | awk -v OFS='\t' '{ print $1"_B", $2}' >> $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.X.txt

module purge
module load R-uoneasy/4.3.2-gfbf-2023a
module load PhyML/3.3.20220408-foss-2023a

## Activate twisst-ete3-p3-6
echo "activate twisst-ete3-p3-6"
conda activate twisst-ete3-p3-6

### Run Genomics general script for calculating trees over a sliding window
## Best to use one thread because it doesnt take that line
## NOTE THIS CODE HAS BEEN MOVED INTO THE MAIN GENOMICS GENERAL DIRECTORY SO IT CAN ACCESS THE GENOMICS.PY SCRIPT

## Run for autosomes
python ~/apps/genomics_general/phyml_sliding_windows.py -T 1 -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
  --prefix $output_dir/${pop1}_${pop2}/${pop1}_${pop2} --windType sites --model GTR --windSize $mywindow -O 0 -M 1 -Ms 1 --indFile $output_dir/${pop1}_${pop2}/ind_file_${pop1}_${pop2}.txt

## Run for PAR
python ~/apps/genomics_general/phyml_sliding_windows.py -T 1 -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz \
  --prefix $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.PAR --windType sites --model GTR --windSize $mywindow -O 0 -M 1 -Ms 1 --indFile $output_dir/${pop1}_${pop2}/ind_file_${pop1}_${pop2}.txt

## Run for X just using female samples
python ~/apps/genomics_general/phyml_sliding_windows.py -T 1 -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz \
  --prefix $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.X --windType sites --model GTR --windSize $mywindow -O 0 -M 1 -Ms 1 --indFile $output_dir/${pop1}_${pop2}/ind_file_${pop1}_${pop2}.X.txt

# CHECK
echo "PHYML RUN SUCCESFULLY"

echo "Print out all populations in input files"
awk '{ print $2 }' $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.txt | sort | uniq > $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}_unique.txt
pop1_sub1=$(awk "FNR==1" $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}_unique.txt)
pop1_sub2=$(awk "FNR==2" $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}_unique.txt)
pop2_sub1=$(awk "FNR==3" $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}_unique.txt)
pop2_sub2=$(awk "FNR==4" $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}_unique.txt)

## Running Twisst (require install of ete3)
python ~/apps/twisst/twisst.py -t $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.trees.gz -w $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.weights.tsv.gz \
       -g ${pop1_sub1} -g ${pop1_sub2} -g ${pop2_sub1} -g ${pop2_sub2} \
       --method ${weight_method} --groupsFile $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.txt

# PAR
python ~/apps/twisst/twisst.py -t $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.PAR.trees.gz -w $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.PAR.weights.tsv.gz \
       -g ${pop1_sub1} -g ${pop1_sub2} -g ${pop2_sub1} -g ${pop2_sub2} \
       --method ${weight_method} --groupsFile $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.txt

# X
python ~/apps/twisst/twisst.py -t $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.X.trees.gz -w $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.X.weights.tsv.gz \
       -g ${pop1_sub1} -g ${pop1_sub2} -g ${pop2_sub1} -g ${pop2_sub2} \
       --method ${weight_method} --groupsFile $output_dir/${pop1}_${pop2}/phased_pop_file_${pop1}_${pop2}.X.txt

conda deactivate
#### Merge all results into one

conda activate genomics-general-p3.13
## Data files
cat $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.data.tsv > $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.data.tsv
awk 'NR != 1 {print}' $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.PAR.data.tsv >> $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.data.tsv
awk 'NR != 1 {print}' $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.X.data.tsv >> $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.data.tsv

## Weights
zcat $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.weights.tsv.gz | bgzip > $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.weights.tsv.gz
zcat $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.PAR.weights.tsv.gz | awk 'NR != 4 {print}' | bgzip >> $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.weights.tsv.gz
zcat $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.X.weights.tsv.gz | awk 'NR != 4 {print}' | bgzip >> $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all.weights.tsv.gz

conda deactivate

## Load R module (DO NOT MOVE TO START OF SCRIPT AS IT BREAKS THE PYTHON VERSION)
module purge
module load R-uoneasy/4.2.1-foss-2022a

## Plotting output of Twisst
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/11.1-twisst.R $output_dir/${pop1}_${pop2}/${pop1}_${pop2}.all

### THEN ON THE CONSOLE (using interactive slurm job) run the code 
### Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/11.3-Twisst-Populations-combined-plot.R 