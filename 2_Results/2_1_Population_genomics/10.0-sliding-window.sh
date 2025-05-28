#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20g
#SBATCH --time=18:00:00
#SBATCH --job-name=sliding-window
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module purge 
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load R-uoneasy/4.2.1-foss-2022a

conda activate genomics-general-p3.13

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
mkdir -p $wkdir/results/sliding-window

## Use Population (waterbody + ecotype)
# grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/species_pairs_sequence_data.csv | 
#    awk -F ',' '{ print $1 " " $10}' > $wkdir/results/sliding-window/pop_file.txt

# Just use Ecotype
grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $13}' > $wkdir/results/sliding-window/pop_file.txt

# Run sliding window script
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
   -o $wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi.csv -f phased -T $SLURM_CPUS_PER_TASK --popsFile $wkdir/results/sliding-window/pop_file.txt -p anad -p resi
## Plot results
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/sliding_window_w25kb_s5kb_m1_Panad_resi"

### Go through each Water body and compare

for waterbody in DUIN LUIB CLAC OBSE
do 
echo "Running comparison just for $waterbody"
# Make folder
mkdir -p $wkdir/results/sliding-window/$waterbody/
## Create pop file bu subsetting
grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
    grep -E "$waterbody" | 
    awk -F ',' '{ print $1, $13}' > $wkdir/results/sliding-window/$waterbody/pop_file_${waterbody}.txt

echo "Caculating Fst"
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
    -o $wkdir/results/sliding-window/$waterbody/sliding_window_w25kb_s5kb_m1_${waterbody}_Panad_resi.csv -f phased -T $SLURM_CPUS_PER_TASK --popsFile $wkdir/results/sliding-window/$waterbody/pop_file_${waterbody}.txt -p anad -p resi

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/$waterbody/sliding_window_w25kb_s5kb_m1_${waterbody}_Panad_resi"

done

## Comparison between CLAC and all other populations
echo "Running comparison just for CLAC"
# Make folder
mkdir -p $wkdir/results/sliding-window/CLAC/

## Create pop file by subsetting (also assign all stream to fw resident)
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
    awk -F ',' '{ print $1, $10}' | sed s/OLST/OLAV/ > $wkdir/results/sliding-window/CLAC/pop_file.txt

for pop in TORM DUIN DUIM LUIM LUIB CLAM OLAV
do
echo "Caculating Fst for CLAC vs ${pop}"
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $wkdir/results/sliding-window/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_${pop}.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/CLAC/pop_file.txt -p CLAC -p ${pop}

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_${pop}"
done

mkdir -p $wkdir/results/sliding-window/OLAV/
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
    awk -F ',' '{ print $1, $10}' | sed s/OLST/OLAV/ > $wkdir/results/sliding-window/OLAV/pop_file.txt

## Comparison between CLAC and other Populations
for pop in OLAM TORM DUIN DUIM LUIM LUIB CLAM
do
echo "Caculating Fst for OLAV vs ${pop}"
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \    -o $wkdir/results/sliding-window/OLAV/sliding_window_w25kb_s5kb_m1_OLAV_${pop}.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/OLAV/pop_file.txt -p OLAV -p ${pop}

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/OLAV/sliding_window_w25kb_s5kb_m1_OLAV_${pop}"
done

############################################
## Comparison between stream and fw resident
############################################

mkdir -p $wkdir/results/sliding-window/stream/
grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
    awk -F ',' '{ print $1, $10}' > $wkdir/results/sliding-window/stream/pop_file.txt

## TORM vs TOST
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_TORM_TOST.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/stream/pop_file.txt -p TORM -p TOST

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_TORM_TOST"

## OLAV vs OLST
python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_OLAV_OLST.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/stream/pop_file.txt -p OLAV -p OLST

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_OLAV_OLST"

# Complete comparison between all stream and fw residents
## OLAV & TORM vs OLST & TOST

grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv |
    awk -F ',' '{ print $1, $10 }' | sed s/TORM/fw/ | sed s/OLAV/fw/ | sed s/OLST/strm/ | sed s/TOST/strm/ > $wkdir/results/sliding-window/stream/pop_file_stream.txt

python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_fw_strm.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/stream/pop_file_stream.txt -p fw -p strm 

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/stream/sliding_window_w25kb_s5kb_m1_fw_strm"

# Complete comparison between all stream and fw residents
## CLAC resi vs all other resi

grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
    awk -F ',' '{ print $1, $10 }' | grep "CLAC" | awk '{print $1, "CLAC_resi" }' > $wkdir/results/sliding-window/CLAC/pop_file.txt

grep -f $wkdir/vcfs/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
    awk -F ',' '{ print $1, $10, $13 }' | grep -E "LUIB|DUIN|OBSE" | awk '{ print $1, $3 }' >> $wkdir/results/sliding-window/CLAC/pop_file.txt

python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
    -o $wkdir/results/sliding-window/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_allresi.csv -f phased -T $SLURM_CPUS_PER_TASK \
    --popsFile $wkdir/results/sliding-window/CLAC/pop_file.txt -p CLAC_resi -p resi 

echo "Plotting"
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.1-sliding-window-plot.R "$wkdir/results/sliding-window/CLAC/sliding_window_w25kb_s5kb_m1_CLAC_resi_allresi"

