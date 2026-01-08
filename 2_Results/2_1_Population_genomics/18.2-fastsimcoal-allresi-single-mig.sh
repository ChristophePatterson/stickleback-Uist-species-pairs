#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH --time=72:00:00
#SBATCH --job-name=fastsimcoal2-allresi-single-mig.sh
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
foldtype=("folded")
# Analyse name
analysis_name=allresi_singlemig_JSFS_Mono_N${SLURM_CPUS_PER_TASK}

## Output
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/${analysis_name}_r${randSNP})
mkdir -p $output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
## Get unique combination of populations
## Get each unique Population (but replacing st with fw)
grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $10, $13 } $13=="st" { print $1, $9, "fw" }' > $output_dir/complete_pop_file.txt

## Subset to samples
awk '$4!="anad" {print $1, $3} $4=="anad" {print $1, "anad"}' $output_dir/complete_pop_file.txt | \
   sort -k 2 > $output_dir/pop_file_unordered.txt
# Get unique waterbodies
# Predefine order of populations
echo -e "CLAC\nLUIB\nOBSE\nDUIN\nanad" > $output_dir/pop_uniq.txt

rm -f $output_dir/pop_file.txt
## loop through each population and set 
for pop in $(cat $output_dir/pop_uniq.txt); do
    echo "Processing population: $pop"
    # Extract individuals for each population
    awk -v popu=$pop '$2==popu {print $1, $2}' $output_dir/pop_file_unordered.txt >> $output_dir/pop_file.txt
done

# Alternative: get unique populations from complete pop file
## awk '{ print $3 }' $output_dir/complete_pop_file.txt | sort | uniq > $output_dir/pop_uniq.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/pop_file.txt > $output_dir/ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites that don't have a at least some (non-zero) minor allele freq, must filter to just snps first.
# With random filtering for reduced input
bcftools view -i 'N_ALT<=1' -S $output_dir/ind_file.txt $vcf | \
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_dir/${analysis_name}_r${randSNP}.vcf.gz

# Copy pre-made vcf 
### cp /gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/allpops_N1_r100000/allpops_N1_r100000.vcf.gz $output_dir/${analysis_name}_r${randSNP}.vcf.gz

## Chosen vcf (used to swap out vcfs in bug testing)
vcf_SFS=$output_dir/${analysis_name}_r${randSNP}

# Calculate number of samples and sequences input to SFS
SAMPcount=$(bcftools query -l $vcf_SFS.vcf.gz | wc -l | awk '{print $1}')
SEQcount=$(bcftools query -l $vcf_SFS.vcf.gz | wc -l | awk '{print $1*2}')

# Calculate number of SNPs input to SFS
SNPcount=$(bcftools view -H $vcf_SFS.vcf.gz | wc -l)

# Deactivate bcftools enviroment
conda deactivate
module purge

######################################
##### Get best proj for easySFS #####
######################################
# Activate easySFS environment
conda activate easySFS-env

echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt --preview > $output_dir/proj_folded.txt

# Display projection preview
cat $output_dir/proj_folded.txt

# Remove old files
rm -f $output_dir/*_proj_long_*.txt
rm -f $output_dir/best_proj_*.txt

## loop through each population and get best projection
for pop in $(cat $output_dir/pop_uniq.txt); do
    echo "Processing population: $pop"
    # Convert projected SFS to long format table
    grep -A 1 $pop $output_dir/proj_$foldtype.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/${pop}_proj_long_$foldtype.txt
    # Extract best projection number
    topprojpop=$(awk '{print $2}' $output_dir/${pop}_proj_long_$foldtype.txt | sort -n | tail -1)
    # Get best proj for each population
    bestprojpop=$(awk -v topproj=$topprojpop '$2==topproj {print $1 }' $output_dir/${pop}_proj_long_$foldtype.txt | sort -n | tail -1)
    echo "Best projection for $pop is $bestprojpop"
    # Save best projection to file
    echo -e "${pop}\t${bestprojpop}" >> $output_dir/best_proj_$foldtype.txt
done

############################
##### Run easySFS #####
############################

## Create folded output
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt -a -f --total-length $SNPcount -o $output_dir/SFS_$foldtype/ --prefix ${analysis_name}_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_dir/best_proj_$foldtype.txt)   

# deactivate easySFS environment
conda deactivate
module purge

##### 
mkdir -p $output_dir/fsc_run
cd $output_dir/fsc_run

## Remove old obs
rm -f ./*.obs
## Create new SFS that has zero in monomorphic sites
# Getting first three lines
awk 'NR<3 {print $0}' $output_dir/SFS_$foldtype/fastsimcoal2/${analysis_name}_${foldtype}_MSFS.obs > ${analysis_name}_${foldtype}_MSFS.obs 
# Cutting first SFS value (monomorphic sites)
echo "$(awk 'NR>=3 {print $0}' $output_dir/SFS_$foldtype/fastsimcoal2/${analysis_name}_${foldtype}_MSFS.obs | cut -d ' ' -f 2-)" >> ${analysis_name}_${foldtype}_MSFS.obs  

# Sum up all values in SFS (non including monomorphic sites)
# Transforming into awk to sum (%.13 increases decimal places)
awk 'NR>=3 {print $0}' ${analysis_name}_${foldtype}_MSFS.obs   | \
    sed 's/ /\n/g' > ${analysis_name}_${foldtype}_MSFS_long.obs
projSFSsites=$(awk -v OFMT=%.13g '{sum += $1} END {print sum}' ${analysis_name}_${foldtype}_MSFS_long.obs)

# Copy joint MAF SFS files 
cp $output_dir/SFS_$foldtype/fastsimcoal2/${analysis_name}_${foldtype}_jointMAF*.obs ./

############################
 ## Create model parameters file ##
############################

echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "5 samples to simulate :" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print $1"$"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//Samples sizes and samples age" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
awk '{print $2}' $output_dir/best_proj_$foldtype.txt >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
awk '{print "0"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
# Historical events
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "4 historical event" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
# Get population numbers for events (remember fsc starts counting at 0)
#OBSE merges into DUIN
echo "TDivResiEast@ $(awk '$1=="OBSE" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE2 0 0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
# CLAC merges into LUIB
echo "TDivResiWest@ $(awk '$1=="CLAC" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="LUIB" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE1 0 0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
# Resi west and east merge
echo "TDivResi@ $(awk '$1=="LUIB" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE3 0 0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs@ $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="anad" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE4 0 0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl

echo "//Number of independent loci [chromosome]" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "1 0" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "1" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl
echo "//per Block:data typ" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl

###### TO FIX ######
echo "FREQ $projSFSsites 0 5.11e-9 OUTEXP"  >> $output_dir/fsc_run/${analysis_name}_${foldtype}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=100
maxNPOP=10000000

echo "// Priors and rules file" > $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "// *********************" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "[PARAMETERS]" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "//all N are in number of haploid individuals" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
# Ecotype ancestral sizes
echo "1 Ancs$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "1 Resi$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
# East and west ancestral sizes
echo "1 ResiWest$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "1 ResiEast$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est

## All popuulations
awk -v mxNPOP=$maxNPOP '{print "1 "$1"$ unif 1000 "mxNPOP" output"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est

# Divergence times
# The lower range limit is an absolute minimum, whereas the upper range is only used as a
# maximum for choosing a random initial value for this parameter. There is actually no upper
# limit to the search range, as this limit can grow by 30% after each cycle
echo "1 TDivAncs@ unif 1000 200000 output" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "1 TDivResi@ unif 1000 TDivAncs@ output paramInRange" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "1 TDivResiWest@ unif 1000 TDivResi@ output paramInRange">> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "1 TDivResiEast@ unif 1000 TDivResi@ output paramInRange">> $output_dir/fsc_run/${analysis_name}_${foldtype}.est

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "0 RESIZE1 = ResiWest$/LUIB$ hide" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "0 RESIZE2 = ResiWest$/DUIN$ hide" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "0 RESIZE3 = Resi$/ResiWest$ hide" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est
echo "0 RESIZE4 = Ancs$/anad$ hide" >> $output_dir/fsc_run/${analysis_name}_${foldtype}.est

## Run fsc
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}_${foldtype}.tpl -n 100000 -e ${analysis_name}_${foldtype}.est -y 4 -m -M -L 40 -c $SLURM_CPUS_PER_TASK > $output_dir/fsc_run/fsc_log_jobID${SLURM_ARRAY_TASK_ID}.txt