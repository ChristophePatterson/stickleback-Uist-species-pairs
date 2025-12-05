#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60g
#SBATCH --time=18:00:00
#SBATCH --array=1-10
#SBATCH --job-name=stairway_plot
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

############################
   # PREPARE ENVIRONMENT #
############################

module purge
source /gpfs01/home/${USER}/.bashrc

# Load easySFS
conda activate easySFS-env

# set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
species=stickleback
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/ploidy_aware_HWEPops_MQ10_BQ20)

## Output
output_dir=($wkdir/results/$vcf_ver/Stairway)
mkdir -p $output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
# echo -e "CLAC\nCLAM\nOBSE\nOBSM\nDUIN\nDUIM\nLUIB\nLUIM\nOLAV\nTORM" > ${output_dir}/pop_list.txt
## Get population that equals slurm array
pop=$(awk -v slurmA=$SLURM_ARRAY_TASK_ID 'NR==slurmA {print $0}' ${output_dir}/pop_list.txt)

## Make directory
mkdir -p $output_dir/$pop

## Use Population (waterbody + ecotype)
grep -w -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_header_2025-04-28.csv | 
   awk -F ',' -v OFS='\t' -v pop=$pop 'NR!=1 && $10==pop { print $1, $10 }' > $output_dir/$pop/${pop}_pop_file.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/$pop/${pop}_pop_file.txt > $output_dir/$pop/${pop}_ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
## With random filtering for reduced input
bcftools view -S $output_dir/$pop/${pop}_ind_file.txt $vcf | \
   bcftools +prune -n 1 -N rand -w 1000bp -O z -o $output_dir/$pop/${pop}_r1000.vcf.gz

# Including all sites
# bcftools view -S $output_dir/$pop/${pop}_ind_file.txt $vcf -O z -o $output_dir/$pop/${pop}.vcf.gz 

## Chosen vcf
vcf_ver=$output_dir/$pop/${pop}_r1000

SAMPcount=$(bcftools query -l $vcf_ver.vcf.gz | wc -l | awk '{print $1}')
SEQcount=$(bcftools query -l $vcf_ver.vcf.gz | wc -l | awk '{print $1*2}')

# Calculate number of SNPs input to SFS
SNPcount=$(bcftools view -H $vcf_ver.vcf.gz | wc -l)

# Deactivate bcftools enviroment
conda deactivate

conda activate easySFS-env
python ~/apps/easySFS/easySFS.py -a -i $vcf_ver.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj_unfolded.txt

# Convert to table
tail -n 2 $output_dir/$pop/${pop}_proj_unfolded.txt | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop}_proj_long_unfolded.txt

# Extract best projection number
topproj=$(awk '{print $2}' $output_dir/$pop/${pop}_proj_long_unfolded.txt | sort -n | tail -1)
bestproj=$(awk -v topproj=$topproj '$2==topproj {print $1 }' $output_dir/$pop/${pop}_proj_long_unfolded.txt | sort -n | tail -1)

# Run easySFS
python ~/apps/easySFS/easySFS.py -i $vcf_ver.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt -a -f --total-length $SNPcount -o $output_dir/$pop/SFS_unfolded/ --prefix ${pop}_unfolded --proj $bestproj

## Input SFS
# get SFS and remove first monomorphic sites
awk 'NR == 3 {print $0}' $output_dir/$pop/SFS_unfolded/fastsimcoal2/${pop}_unfolded_MSFS.obs | cut -d ' ' -f 2-$(expr $bestproj) > $output_dir/$pop/${pop}_input_sfs.obs

# Deactivate easySFS
conda deactivate

## Purge modules
module purge

# Load java
module load java-uoneasy/17.0.6

## Stairway filepath
stairpath=/gpfs01/home/mbzcp2/apps/stairway-plot-v2/stairway_plot_v2.2

# Remove prior results
rm -r $stairpath/$pop/

############################
## Create blue print file ##
############################

echo "#Blueprint file for $pop" > $stairpath/$pop.blueprint.txt
echo "#input setting" >> $stairpath/$pop.blueprint.txt
echo "popid: $pop"  >> $stairpath/$pop.blueprint.txt # id of the population (no white space)"
# echo "nseq: $(wc -l $stairpath/${pop}_ind_file.txt| awk '{print $1}')" >> $stairpath/$pop.blueprint.txt # number of sequences
echo "nseq: $bestproj" >> $stairpath/$pop.blueprint.txt # number of sequences
echo "L: $SNPcount" >> $stairpath/$pop.blueprint.txt # total number of observed nucleic sites, including polymorphic and monomorphic
echo "whether_folded: false" >> $stairpath/$pop.blueprint.txt # whethr the SFS is folded (true or false)
echo "SFS: $(awk 'NR == 1 {print $0}' $output_dir/$pop/${pop}_input_sfs.obs)" >> $stairpath/$pop.blueprint.txt # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
# Parameters
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 29 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
echo "pct_training: 0.67" >> $stairpath/$pop.blueprint.txt # percentage of sites for training
echo " nrand: 7 15 22 28" >> $stairpath/$pop.blueprint.txt # number of random break points for each try (separated by white space)
echo "project_dir: $pop" >> $stairpath/$pop.blueprint.txt # project directory
echo "stairway_plot_dir: $stairpath/stairway_plot_es" >> $stairpath/$pop.blueprint.txt # directory to the stairway plot files
echo "ninput: 200" >> $stairpath/$pop.blueprint.txt # number of input files to be created for each estimation
echo "theta_upper_bound: 0.2" >> $stairpath/$pop.blueprint.txt # the maximum of theta used in the search algorithm, default is 0.2. Increase this value only if seeing confident intervals converged to the same value as a plateau on the plot.
echo "dimension_factor: 2000" >> $stairpath/$pop.blueprint.txt # this parameter determin the maximum number of iteration of the search algorithm, default is 2000.
echo "#random_seed: 6" >> $stairpath/$pop.blueprint.txt
echo "#output setting" >> $stairpath/$pop.blueprint.txt
echo "mu: 5.11e-9" >> $stairpath/$pop.blueprint.txt # assumed mutation rate per site per generation" 
echo "year_per_generation: 1" >> $stairpath/$pop.blueprint.txt # assumed generation time (in years)
echo "#plot setting" >> $stairpath/$pop.blueprint.txt
echo "plot_title: $pop" >> $stairpath/$pop.blueprint.txt # title of the plot
echo "xrange: 0.1,10000" >> $stairpath/$pop.blueprint.txt # Time (1k year) range; format: xmin,xmax; "0,0" for default
echo "yrange: 0,0" >> $stairpath/$pop.blueprint.txt # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
echo "xspacing: 2" >> $stairpath/$pop.blueprint.txt # X axis spacing
echo "yspacing: 2" >> $stairpath/$pop.blueprint.txt # Y axis spacing
echo "fontsize: 12" >> $stairpath/$pop.blueprint.txt # Font size

cd $stairpath
## Run Stairway plot
java -cp stairway_plot_es Stairbuilder $stairpath/$pop.blueprint.txt

bash $stairpath/$pop.blueprint.txt.sh

## Copy results into results folder
mkdir -p $output_dir/results_unfolded/
cp $stairpath/$pop/${pop}* $output_dir/results_unfolded/

