#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1g
#SBATCH --array=1-66
#SBATCH --time=04:00:00
#SBATCH --job-name=sliding-window-pops
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## Due to overlap in writing files, if slurm array is not equal to 1 then wait 30 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

############################
   # PREPARE ENVIRONMENT #
############################
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
vcf_ver=ploidy_aware_HWEPops

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/$vcf_ver/sliding-window/

# Make folder
mkdir -p $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison

## Create popfile and unique combination of population (if doesn't exist already)
if [ ! -f $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt ]; then
    grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
        awk -F ',' '{ print $1, $10}' | sed s/OLST/OLAV/ > $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt
        ## Add outgroup for CROC
    grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
        awk -F ',' '$10=="CROC" || $10=="TORM" {print $1, $10} $10=="OLAV" || $10=="OLST" {print $1, "OLAV"} $10=="OLAM" {print $1, $10}' >> $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt

    # Get all unique populations
    awk '{print $2}' $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt | sort | uniq > $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt
    ## Get all unique combination of populations
    Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/Helper_scripts/Get_combinations.R $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt 2
    wc -l $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt 
fi

## Get unique combination of waterbodies
pop1=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $1 }' $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt)
pop2=$(awk -v slurmID="$SLURM_ARRAY_TASK_ID" 'FNR == slurmID { print $2 }' $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file_uniq.txt_combn.txt)

# Create output directory
outputdir=($wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/${pop1}_${pop2})
mkdir -p $outputdir

module purge
conda activate genomics-general-p3.13

echo "Running Fst sliding window for ${pop1} and ${pop2}."

### Run sliding window for all autosomes
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
   -o $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_auto.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt -p ${pop1} -p ${pop2}

# For PAR
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz \
   -o $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_PAR.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt -p ${pop1} -p ${pop2}

# For X
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_SNPs.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz \
   -o $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_X.csv -f phased --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt -T $SLURM_CPUS_PER_TASK \
   --popsFile $wkdir/results/$vcf_ver/sliding-window/All_Pop_comparison/pop_file.txt -p ${pop1} -p ${pop2}

# Merge all Fst calculations into a single file
cat $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_auto.csv > $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}.csv
awk FNR!=1 $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_PAR.csv >> $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}.csv
awk FNR!=1 $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}_X.csv >> $outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}.csv

##### ## Plot results
##### ## Deactivate conda and load R
conda deactivate
module load R-uoneasy/4.2.1-foss-2022a

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$outputdir/sliding_window_w25kb_s5kb_m1_${pop1}_${pop2}"
