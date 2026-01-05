#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60g
#SBATCH --time=24:00:00
#SBATCH --job-name=fastsimcoal2-allpops
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
randSNP=100000

# folded or unfold
foldtype=("folded")
# Analyse name
analysis_name=allpops

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
# Get unique waterbodies
awk '{ print $3 }' $output_dir/complete_pop_file.txt | sort | uniq > $output_dir/waterbody_uniq.txt

## Subset to samples
awk '{print $1, $3}' $output_dir/complete_pop_file.txt | \
   sort -k 2 > $output_dir/pop_file.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/complete_pop_file.txt > $output_dir/ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites that don't have a at least some (non-zero) minor allele freq, must filter to just snps first.
# With random filtering for reduced input
bcftools view -v snps -i 'N_ALT=1' -S $output_dir/ind_file.txt $vcf | \
    bcftools +fill-tags -- -t AN,AC,AF,MAF | \
    bcftools view -q '0.00000001:minor' -Q '0.9999999:minor' |
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_dir/${analysis_name}_r${randSNP}.vcf.gz

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
for pop in $(cat $output_dir/waterbody_uniq.txt); do
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

pop0=$(awk 'NR==1 {print $2}' $output_dir/best_proj_$foldtype.txt)
############################
##### Run easySFS #####
############################

## Create folded output
if [[ $foldtype == "folded" ]]; then
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt -a -f --total-length $SNPcount -o $output_dir/SFS_$foldtype/ --prefix ${analysis_name}_$foldtype \
--proj=$(awk 'NR==1 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==6 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==7 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==8 {print $2}' $output_dir/best_proj_$foldtype.txt)
   
fi

##### 
mkdir -p $output_dir/fsc_run
cd $output_dir/fsc_run

# Sum of all polymorphic sites included in SFS
# Getting third line
# Cutting first SFS value (monomorphic sites)
# Transforming into awk to sum (%.13 increases decimal places)

## Remove old obs
rm -f ./*.obs
## Create new SFS that has zero in monomorphic sites
awk 'NR<3 {print $0}' $output_dir/$pop/SFS_$foldtype/fastsimcoal2/${pop}_${foldtype}_MSFS.obs > ${pop0}-${pop1}_MSFS.obs 
echo "0 $(awk 'NR==3 {print $0}' $output_dir/$pop/SFS_$foldtype/fastsimcoal2/${pop}_${foldtype}_MSFS.obs | cut -d ' ' -f 2-)" >> ${pop0}-${pop1}_MSFS.obs 

# Sum up all values in SFS (non including monomorphic sites)
projSFSsites=$(awk 'NR==3 {print $0}' ${pop0}-${pop1}_MSFS.obs | \
    sed 's/ /\n/g' | \
    awk -v OFMT=%.13g '{sum += $1} END {print sum}')

cat $output_dir/waterbody_uniq.txt
## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/fsc_run/${analysis_name}.tpl
echo "8 samples to simulate :" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/fsc_run/${analysis_name}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print "NPOP"$1}' $output_dir/waterbody_uniq.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Samples sizes and samples age" >> $output_dir/fsc_run/${analysis_name}.tpl
awk '{print $2}' $output_dir/best_proj_$foldtype.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/fsc_run/${analysis_name}.tpl
awk '{print "0"}' $output_dir/waterbody_uniq.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "0" >> $output_dir/fsc_run/${analysis_name}.tpl

$(awk '$1=="LUIB" {print NR}' $output_dir/waterbody_uniq.txt)

echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "7 historical event" >> $output_dir/fsc_run/${analysis_name}.tpl
# CLAC merges into LUIB
echo "TDivRWest $(awk '$1=="CLAC" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="LUIB" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE1 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
#OBSE merges into DUIN
echo "TDivREast $(awk '$1=="OBSE" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="DUIN" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE2 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
# CLAM merges into LUIM
echo "TDivMWest $(awk '$1=="CLAM" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="LUIM" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE3 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
# OBSM merges into DUIM
echo "TDivMEast $(awk '$1=="OBSM" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="DUIM" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE4 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
# Resi west and east merge
echo "TDivResi $(awk '$1=="LUIB" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="DUIN" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE5 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
# Migration west and east merge
echo "TDivMigr $(awk '$1=="LUIM" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="DUIM" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE6 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs $(awk '$1=="DUIN" {print NR}' $output_dir/waterbody_uniq.txt) $(awk '$1=="DUIM" {print NR}' $output_dir/waterbody_uniq.txt) 1 RESIZE7 0 1" >> $output_dir/fsc_run/${analysis_name}.tpl

echo "//Number of independent loci [chromosome]" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 0" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//per Block:data typ" >> $output_dir/fsc_run/${analysis_name}.tpl

###### TO FIX ######
echo "FREQ $projSFSsites 0 5.11e-9 OUTEXP"  >> $output_dir/fsc_run/${analysis_name}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=1000
maxNPOP=1000000

echo "// Priors and rules file" > $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "// *********************" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "[PARAMETERS]" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "//all N are in number of haploid individuals" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
## All popuulations
awk -v mxNPOP=$maxNPOP '{print "1 NPOP"$1" unif 1000 "mxNPOP" output"}' $output_dir/waterbody_uniq.txt >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
# East and west ancestral sizes
echo "1 ResiWest unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 ResiEast unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 MigrWest unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 MigrEast unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
# Ecotype ancestral sizes
echo "1 Resi 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 Migr unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 Ancs unif 1000 $maxNPOP output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est

# Divergence times
echo "1 TDivRWest unif 1000 TDivResi output paramInRange">> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivREast unif 1000 TDivResi output paramInRange">> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivMWest unif 1000 TDivMigr output paramInRange">> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivMEast unif 1000 TDivMigr output paramInRange">> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivResi unif 1000 TDivAncs output paramInRange" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivMigr unif 1000 TDivAncs output paramInRange" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 TDivAncs unif 1000 200000 output" >> $output_dir/fsc_run/${analysis_name}.tpl

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE1 = ResiWest/NPOPLUIB hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE2 = ResiEast/NPOPDUIN hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE3 = MigrWest/NPOPLUIM hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE4 = MigrEast/NPOPDUIM hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE5 = Resi/ResiEast hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE6 = Migr/MigrEast hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE7 = Ancs/Migr hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est


## Run fsc
# ~/apps/fsc28_linux64/fsc28 -t $pop0-$pop1.tpl -n 1000000 -e $pop0-$pop1.est -m -u -M -L 1000 -c $SLURM_CPUS_PER_TASK -q 