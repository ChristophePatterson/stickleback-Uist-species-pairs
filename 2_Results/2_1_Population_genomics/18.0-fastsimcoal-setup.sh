#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=10:00:00
#SBATCH --array=1-100
#SBATCH --job-name=fastsimcoal2-multi-model-comparison
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)
randSNP=10000

# folded or unfold
foldtype=("unfolded")
# Global output directory
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/model_selection_${foldtype}_r${randSNP})
SFS_name=(SFS_$SLURM_ARRAY_TASK_ID)
output_SFS_dir=($output_dir/SFS/${SFS_name})

# Create output directories
mkdir -p $output_dir
mkdir -p $output_SFS_dir

####################################
   # GENERATE RANDOM SNP DATASETS #
####################################

## Input master vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
## Get unique combination of populations
## Get each unique Population (but replacing st with fw)
grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $10, $13 } $13=="st" { print $1, $9, "fw" }' > $output_SFS_dir/complete_pop_file.txt

# Create file with list of individuals
awk '{print $1}' $output_SFS_dir/complete_pop_file.txt > $output_SFS_dir/ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites where the altnative allele is fixed in all populations, which is illegal for fastsimcoal2 interpretation of alleles in coalescent theory
# With random filtering for reduced input
bcftools view -i 'N_ALT<=1' -S $output_SFS_dir/ind_file.txt $vcf | \
    bcftools +fill-tags -- -t AN,AC,AF,MAF | \
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_SFS_dir/${SFS_name}.vcf.gz

# If foldtype is unfolded filter out all alt fixed sites
if [[ $foldtype == "unfolded" ]]; then
    bcftools view -e 'N_ALT=1 && AC=AN' $output_SFS_dir/${SFS_name}.vcf.gz -O z -o $output_SFS_dir/${SFS_name}_filtNAlt1.vcf.gz
    # Overwrite original vcf with filtered vcf
    mv $output_SFS_dir/${SFS_name}_filtNAlt1.vcf.gz $output_SFS_dir/${SFS_name}.vcf.gz
    # Remove intermediate file
    rm -f $output_SFS_dir/${SFS_name}_filtNAlt1.vcf.gz
fi

########################################
   # Generate SFS for all populations #
########################################

# Calculate number of samples and sequences input to SFS
SAMPcount=$(bcftools query -l $output_SFS_dir/${SFS_name}.vcf.gz | wc -l | awk '{print $1}')
SEQcount=$(bcftools query -l $output_SFS_dir/${SFS_name}.vcf.gz | wc -l | awk '{print $1*2}')

# Calculate number of SNPs input to SFS
SNPcount=$(bcftools view -H $output_SFS_dir/${SFS_name}.vcf.gz | wc -l)

# Deactivate bcftools enviroment
conda deactivate
module purge

##### Get best proj for easySFS #####
# Activate easySFS environment
conda activate easySFS-env

# Predefine order of populations
echo -e "CLAC\nLUIB\nOBSE\nDUIN\nCLAM\nLUIM\nOBSM\nDUIM" > $output_SFS_dir/pop_all_uniq.txt

rm -f $output_SFS_dir/pop_all_file.txt
## loop through each population and set 
for pop in $(cat $output_SFS_dir/pop_all_uniq.txt); do
    echo "Processing population: $pop"
    # Extract individuals for each population
    awk -v popu=$pop '$3==popu {print $1, $3}' $output_SFS_dir/complete_pop_file.txt >> $output_SFS_dir/pop_all_file.txt
done

# Create projection preview file
# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $output_SFS_dir/${SFS_name}.vcf.gz --unfolded -p $output_SFS_dir/pop_all_file.txt --preview > $output_SFS_dir/${SFS_name}_all_proj_$foldtype.txt
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_all_file.txt --preview > $output_SFS_dir/${SFS_name}_all_proj_$foldtype.txt
fi

# Display projection preview
cat $output_SFS_dir/${SFS_name}_all_proj_$foldtype.txt

# Remove old files
rm -f $output_SFS_dir/${SFS_name}_all_proj_long_$foldtype.txt
rm -f $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt

## loop through each population and get best projection
for pop in $(cat $output_SFS_dir/pop_all_uniq.txt); do
    echo "Processing population: $pop"
    # Convert projected SFS to long format table
    grep -A 1 $pop $output_SFS_dir/${SFS_name}_all_proj_$foldtype.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_SFS_dir/${SFS_name}_${pop}_all_proj_long_$foldtype.txt
    # Extract best projection number
    topprojpop=$(awk '{print $2}' $output_SFS_dir/${SFS_name}_${pop}_all_proj_long_$foldtype.txt | sort -n | tail -1)
    # Get best proj for each population
    bestprojpop=$(awk -v topproj=$topprojpop '$2==topproj {print $1 }' $output_SFS_dir/${SFS_name}_${pop}_all_proj_long_$foldtype.txt | sort -n | tail -1)
    echo "Best projection for $pop is $bestprojpop"
    # Save best projection to file
    echo -e "${pop}\t${bestprojpop}" >> $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt
done

grep -A 1 /gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/model_selection_unfolded_r10000/SFS/SFS_1/SFS_1_all_proj_unfolded.txt


## Create folded output
# Create SFS
# runs code but cancels if projection takes longer than 30 seconds. 
# jointMAF are produced quickly but MSFS files can take a long time to produce for large datasets, which are not needed here.
# There is no hope of making a multiSFS with this number of populations, it would be huge (~32 million cells)

# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
timeout 120s \
python ~/apps/easySFS/easySFS.py -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_all_file.txt --unfolded -v -a -f --total-length $SNPcount -o $output_SFS_dir/SFS_all_$foldtype/ --prefix ${SFS_name}_all_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==6 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==7 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==8 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt)
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
timeout 120s \
python ~/apps/easySFS/easySFS.py -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_all_file.txt -v -a -f --total-length $SNPcount -o $output_SFS_dir/SFS_all_$foldtype/ --prefix ${SFS_name}_all_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==6 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),\
$(awk 'NR==7 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt),$(awk 'NR==8 {print $2}' $output_SFS_dir/${SFS_name}_all_best_proj_$foldtype.txt)
fi

###################################################
   # Generate SFS for all resi but single anad #
###################################################

##### Get best proj for easySFS #####
# Activate easySFS environment

## Predefine populations but with anad as separate population
awk '$4!="anad" {print $1, $3} $4=="anad" {print $1, "anad"}' $output_SFS_dir/complete_pop_file.txt | \
   sort -k 2 > $output_SFS_dir/pop_file_unordered.txt

# Predefine order of populations
echo -e "CLAC\nLUIB\nOBSE\nDUIN\nanad" > $output_SFS_dir/pop_sigMig_uniq.txt

rm -f $output_SFS_dir/pop_sigMig_file.txt

## loop through each population and set 
for pop in $(cat $output_SFS_dir/pop_sigMig_uniq.txt); do
    echo "Processing population: $pop"
    # Extract individuals for each population
    awk -v popu=$pop '$2==popu {print $1, $2}' $output_SFS_dir/pop_file_unordered.txt >> $output_SFS_dir/pop_sigMig_file.txt
done


# Create projection preview file
# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $output_SFS_dir/${SFS_name}.vcf.gz --unfolded -p $output_SFS_dir/pop_sigMig_file.txt --preview > $output_SFS_dir/${SFS_name}_sigMig_proj_$foldtype.txt
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_sigMig_file.txt --preview > $output_SFS_dir/${SFS_name}_sigMig_proj_$foldtype.txt
fi

# Display projection preview
cat $output_SFS_dir/${SFS_name}_sigMig_proj_$foldtype.txt

# Remove old files
rm -f $output_SFS_dir/${SFS_name}_*_sigMig_proj_long_$foldtype.txt
rm -f $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt

## loop through each population and get best projection
for pop in $(cat $output_SFS_dir/pop_sigMig_uniq.txt); do
    echo "Processing population: $pop"
    # Convert projected SFS to long format table
    grep -A 1 $pop $output_SFS_dir/${SFS_name}_sigMig_proj_$foldtype.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_SFS_dir/${SFS_name}_${pop}_sigMig_proj_long_$foldtype.txt
    # Extract best projection number
    topprojpop=$(awk '{print $2}' $output_SFS_dir/${SFS_name}_${pop}_sigMig_proj_long_$foldtype.txt | sort -n | tail -1)
    # Get best proj for each population
    bestprojpop=$(awk -v topproj=$topprojpop '$2==topproj {print $1 }' $output_SFS_dir/${SFS_name}_${pop}_sigMig_proj_long_$foldtype.txt | sort -n | tail -1)
    echo "Best projection for $pop is $bestprojpop"
    # Save best projection to file
    echo -e "${pop}\t${bestprojpop}" >> $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt
done

## Create folded output
# Create SFS
# runs code but cancels if projection takes longer than 30 seconds. 
# jointMAF are produced quickly but MSFS files can take a long time to produce for large datasets, which are not needed here.
# There is no hope of making a multiSFS with this number of populations, it would be huge (~32 million cells)

# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
timeout 900s \
python ~/apps/easySFS/easySFS.py -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_sigMig_file.txt --unfolded -v -a -f --total-length $SNPcount -o $output_SFS_dir/SFS_sigMig_$foldtype/ --prefix ${SFS_name}_sigMig_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt)
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
timeout 900s \
python ~/apps/easySFS/easySFS.py -i $output_SFS_dir/${SFS_name}.vcf.gz -p $output_SFS_dir/pop_sigMig_file.txt -v -a -f --total-length $SNPcount -o $output_SFS_dir/SFS_$foldtype/ --prefix ${SFS_name}_sigMig_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_SFS_dir/${SFS_name}_sigMig_best_proj_long_$foldtype.txt)
fi



