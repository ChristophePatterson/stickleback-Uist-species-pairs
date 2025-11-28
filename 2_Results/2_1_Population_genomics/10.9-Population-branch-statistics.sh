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
#SBATCH --job-name=Population-branch-statistics
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################
source /gpfs01/home/${USER}/.bashrc
# load modules
module purge 
conda activate genomics-general-p3.13

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

# Using scripts from https://github.com/simonhmartin/genomics_general?tab=readme-ov-file

## Create input pop file
output_dir=($wkdir/results/$vcf_ver/sliding-window/PopBranchStats)
mkdir -p ${output_dir}

## Window size in Kb
wndsize=(100)
sldsize=(50)

## Use Population (waterbody + ecotype)
grep -w -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv | 
   awk -F ',' -v OFS='\t' 'NR!=1 && $13!="st" { print $1, $13 } NR!=1 && $13=="st" { print $1, "fw" }' > ${output_dir}/pop_file.txt

# Just use Ecotype
# grep -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
#    awk -F ',' -v OFS='\t' '{ print $1, $13}' > ${output_dir}/pop_file.txt

# Female file just for X
grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt ${output_dir}/pop_file.txt > ${output_dir}/pop_file_females.txt
# Pop file just for Y
grep -w -f $wkdir/vcfs/$vcf_ver/male_samples.txt ${output_dir}/pop_file.txt > ${output_dir}/pop_file_males.txt

# Run PopDist window script
# For autosomessldzize
python ~/apps/genomics_general/popgenWindows.py -w ${wndsize}000 -s ${sldsize}000 -m 1 --analysis popDist popPairDist indHet -g $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.geno.gz \
   -o ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_auto.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile ${output_dir}/pop_file.txt -p resi -p anad -p fw

# For PAR
python ~/apps/genomics_general/popgenWindows.py -w ${wndsize}000 -s ${sldsize}000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP5.MEANGTDP5_200.Q60.MAF2.PAR.geno.gz \
   -o ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_PAR.csv -f phased --ploidy 2 -T $SLURM_CPUS_PER_TASK \
   --popsFile ${output_dir}/pop_file.txt -p resi -p anad -p fw

# For X 
python ~/apps/genomics_general/popgenWindows.py -w ${wndsize}000 -s ${sldsize}000 -m 1 --analysis popDist popPairDist -g $wkdir/vcfs/$vcf_ver/${species}_all.NOGTDP2.MEANGTDP2_200.Q60.MAF2.X.geno.gz \
   -o ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_X.csv -f phased --ploidyFile $wkdir/vcfs/$vcf_ver/ploidy_X.txt -T $SLURM_CPUS_PER_TASK \
   --popsFile ${output_dir}/pop_file_females.txt -p resi -p anad -p fw

## Remove individual het from autosome calculations so can merge with X samples
cut -d, -f$(
  head -1 ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_auto.csv | tr ',' '\n' | nl -w1 -s: |
  grep -v "het_" | cut -d: -f1 | tr '\n' ',' | sed 's/,$//'
) ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_auto.csv > ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_auto_pop.csv

## Merge files together
awk -F ',' -v OFS=',' '{print $0,"MF"}' ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_auto_pop.csv > ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_APARX.csv
awk -F ',' -v OFS=',' 'FNR!=1 {print $0,"MF"}' ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_PAR.csv >> ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_APARX.csv
awk -F ',' -v OFS=',' 'FNR!=1 {print $0,"F"}' ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_X.csv >> ${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair_APARX.csv

## Deactivate conda
conda deactivate
# Load R
module load R-uoneasy/4.2.1-foss-2022a
## Plot results
Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.7-sliding-window-PopHet-plot.R "${output_dir}/sliding_window_w${wndsize}kb_s${sldsize}kb_m1_PopPair"

