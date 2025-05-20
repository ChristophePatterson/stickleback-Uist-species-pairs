#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH --time=18:00:00
#SBATCH --job-name=SNP-PCA-LEA
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out
  
############################
   # PREPARE ENVIRONMENT #
############################

# load modules
module load R-uoneasy/4.2.1-foss-2022a
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/SambaR
mkdir -p $outdir

## Plink/Sambar cant have samples with "_" - replace with "-"
bcftools query -l $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz > $outdir/samples.txt
sed s/_/-/ $outdir/samples.txt > $outdir/samples_recode.txt
# Reheader samples
bcftools reheader -s $outdir/samples_recode.txt $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.vcf.gz > $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.vcf.gz

## Convert to Plink format to inlcude input into SambaR
mkdir -p $wkdir/vcfs/plink
plink -vcf  $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader
plink --file $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/plink/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair-wOG.reheader.recode

## Create input popfile for SambaR
echo -e "name\tpop" > $outdir/pop_file_Population.txt
grep -f $outdir/samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $1, $10}' | sed s/NA/Lubec/ | sed s/_/-/ >> $outdir/pop_file_Population.txt

## Run SambaR
# Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/09.0-SambaR.R

################################
  # Genomic sex determination  # 
################################

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/Sex_determination.R

########################
  # LEA - PCA & SNMF  # 
########################

## Custom analysis

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/09-SNP-analysis.R

#####################
  # fastStructure  # 
#####################

## conda create -n faststucture-p2 python=2
module purge
conda activate faststucture-p2
module load bcftools-uoneasy/1.18-GCC-13.2.0
module load plink-uoneasy/2.00a3.7-foss-2023a-highcontig

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/faststructure
mkdir -p $outdir

## Create input file for fasstructure
inputfile=(stickleback_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.rand1000)

## Remove non paired species 
bcftools view -S $wkdir/vcfs/${species}_subset_samples.txt $wkdir/vcfs/${inputfile}.vcf.gz | \
    bcftools view --min-ac 2[minor] -O z -o $wkdir/vcfs/${inputfile}.rand1000.SpPair.vcf.gz

## Convert to plink format
plink -vcf  $wkdir/vcfs/${inputfile}.rand1000.SpPair.vcf.gz --allow-extra-chr -recode --out $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair
plink --file $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair --chr-set 95 --allow-extra-chr --make-bed --recode A --out $wkdir/vcfs/plink/${inputfile}.rand1000.SpPair.recode

## Create popfile for plotting results
grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '{ print $10 }' > $outdir/popfile.txt

## Loop through running faststructure with K= 1 to 5
# Set i as 0
i=0
## Start looop
while [ $i -ne 20 ]
do
        ## Add 1 to K
        i=$(($i+1))
        echo "Running faststructre for K-$i"
        ## Run structure
        structure.py --input=$wkdir/vcfs/plink/${inputfile}.rand1000.SpPair.recode \
            --format bed -K $i --output=$outdir/${species}_SpPair \
            --cv=20 --full

        ## Plot stucture make sure you have added code to distruct.py
        # beforethe to the code "import matplotlib.pyplot as plot"
        # #import matplotlib as mpl
        # #mpl.use('Agg')

        distruct.py --input=$outdir/${species}_SpPair \
             -K $i --output=$outdir/${species}_SpPair_K${i}.pdf --popfile=$outdir/popfile.txt \
             --title=${species}_SpPair_K${i}

done

conda deactivate

## Sliding window FST using genomics general

conda activate genomics-general-p3.13

outdir=/gpfs01/home/mbzcp2/data/sticklebacks/results/popgen
mkdir -p $outdir

## Create popfile
## grep -f $wkdir/vcfs/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
##     awk -F ',' -v OFS='\t' '{ print $1, $10}' > $outdir/pop_file_Population.txt
## 
## python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis indHet -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
##    -o $outdir/sliding_window_w25kb_s5kb_m1_Popgen.csv -f phased -T $SLURM_CPUS_PER_TASK \
##    --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 
## 
## python ~/apps/genomics_general/popgenWindows.py -w 25000 -s 5000 -m 1 --analysis hapStats -g $wkdir/vcfs/${species}_SNPs.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2_SpPair.geno.gz \
##    -o $outdir/sliding_window_w25kb_s5kb_m1_hapStats.csv -f phased -T $SLURM_CPUS_PER_TASK \
##    --popsFile $outdir/pop_file_Population.txt -p CLAC -p CLAM -p DUIM -p DUIN -p LUIB -p LUIM -p OBSE -p OBSM 


conda deactivate