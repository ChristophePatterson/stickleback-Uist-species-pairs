#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=12:00:00
#SBATCH --array=1-21
#SBATCH --job-name=sliding-window-pca
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## Due to overlap in writing files, if slurm array is not equal to 1 then wait 15 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc
conda activate bcftools-env
module load R-uoneasy/4.2.1-foss-2022a

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

########################
  # LEA - PCA & SNMF  # 
########################

## Run variables
vcf=${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.AX
vcf_full=$wkdir/vcfs/$vcf_ver/$vcf.vcf.gz
wndsize=25000
wndslid=5000
run_analysis="TRUE"

output_dir=/gpfs01/home/mbzcp2/data/sticklebacks/results/$vcf_ver/sliding-window/pca/Anad_resi/wndsize${wndsize}_wndslid${wndslid}
mkdir -p $output_dir

## Create config files if this is the first array
if [ $SLURM_ARRAY_TASK_ID == "1" ]; then
   ## Get list of chromosomes to use
   bcftools query -f '%CHROM\n' $vcf_full | sort | uniq > $output_dir/chrom_list.txt 
   ### Populations to use
   bcftools query -l $vcf_full > $output_dir/samples_in_vcf.txt
   echo -e "OBSE\nDUIN\nLUIB\nCLAC" > $output_dir/Pops_interest.txt
   #### Extract sample information
   awk -F ',' -v OFS='\t' '{ print $1, $13 ,$9}' /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | \
         grep -f $output_dir/samples_in_vcf.txt | \
         grep -f $output_dir/Pops_interest.txt > $output_dir/pop_file.txt
   awk '{print $1}' $output_dir/pop_file.txt > $output_dir/samples.txt
fi

# Get chromosome number for array
chr=$(awk "FNR==$SLURM_ARRAY_TASK_ID" $output_dir/chrom_list.txt)
## Create output director
mkdir -p $output_dir/$chr

# Copy over sample file so can be set for each chromosome specifically
cat $output_dir/samples.txt > $output_dir/$chr/samples.txt

## If chr is equal to X chromosome only include female samples
if [ $chr == 'NC_053230.1' ]; then
   echo "$chr is sex chromosome so subsetting just to Female samples"
   grep -w -f $wkdir/vcfs/$vcf_ver/female_samples.txt $output_dir/samples.txt > $output_dir/$chr/samples.txt
fi

# Subset to specific chromosome
bcftools view -r $chr -S $output_dir/$chr/samples.txt --min-ac 2:minor -O z -o $output_dir/$chr/stickleback.$chr.vcf.gz $vcf_full

# Remove result if previously created
if [ $run_analysis == "TRUE" ]; then
   rm -f $output_dir/$chr/stickleback.${chr}_sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt

   # AND Run sliding window PCA
   Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca.R \
        $output_dir/$chr/stickleback.$chr.vcf.gz $vcf_ver $wndsize $wndslid $run_analysis
fi

# Remove subset vcf
rm -f $output_dir/$chr/stickleback.$chr.vcf.gz
# Remove temp PCA files
rm -f $output_dir/$chr/*.lfmm
rm -f $output_dir/$chr/*.geno
rm -f -r $output_dir/$chr/*.pca
rm -f $output_dir/$chr/*.pcaProject

## Merge all PCA and MDS files together - if there are 21 files already created

pcafilesNo=$(ls $output_dir/*/stickleback.*_sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt  | wc -l)
if [ $pcafilesNo == 21 ]; then
   sleep 10
   echo "All $pcafilesNo, perm files created so merging output from all"
   echo -e "sample,chr,start,end,nsnps,nsamps,PCA1,PCA2,MDS1,MDS2" > $output_dir/sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt
   awk FNR!=1 $output_dir/*/stickleback.*_sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt >> $output_dir/sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt
   ## Plot in R
   Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/10.6-sliding-window-pca-plot.R \
      $output_dir/ sliding-window_pca_wndsize${wndsize}_wndslid${wndslid}.txt
else
   echo "There are only $pcafilesNo permutation files so not merging"
fi


