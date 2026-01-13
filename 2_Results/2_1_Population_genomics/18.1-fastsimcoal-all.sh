#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60g
#SBATCH --time=48:00:00
#SBATCH --array=1-100
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
randSNP=10000

# folded or unfold
foldtype=("unfolded")
# Analyse name
analysis_dir=allpops_r${randSNP}_Mono_filtNAlt1_${foldtype}_N12_bootstrap
analysis_name=allpops_r${randSNP}_Mono_filtNAlt1_${foldtype}_N12_A${SLURM_ARRAY_TASK_ID}

## Output
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/$analysis_dir/${analysis_name})
mkdir -p $output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
## Get unique combination of populations
## Get each unique Population (but replacing st with fw)
grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
    awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $10, $13 } $13=="st" { print $1, $9, "fw" }' > $output_dir/complete_pop_file.txt
# Get unique waterbodies
# Predefine order of populations
echo -e "CLAC\nLUIB\nOBSE\nDUIN\nCLAM\nLUIM\nOBSM\nDUIM" > $output_dir/pop_uniq.txt

rm -f $output_dir/pop_file.txt
## loop through each population and set 
for pop in $(cat $output_dir/pop_uniq.txt); do
    echo "Processing population: $pop"
    # Extract individuals for each population
    awk -v popu=$pop '$3==popu {print $1, $3}' $output_dir/complete_pop_file.txt >> $output_dir/pop_file.txt
done

# Alternative: get unique populations from complete pop file
## awk '{ print $3 }' $output_dir/complete_pop_file.txt | sort | uniq > $output_dir/pop_uniq.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/pop_file.txt > $output_dir/ind_file.txt

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites where the altnative allele is fixed in all populations, which is illegal for fastsimcoal2 interpretation of alleles in coalescent theory
# With random filtering for reduced input
bcftools view -i 'N_ALT<=1' -S $output_dir/ind_file.txt $vcf | \
    bcftools +fill-tags -- -t AN,AC,AF,MAF | \
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_dir/${analysis_name}.vcf.gz

# If foldtype is unfolded filter out all alt fixed sites
if [[ $foldtype == "unfolded" ]]; then
    bcftools view -e 'N_ALT=1 && AC=AN' $output_dir/${analysis_name}.vcf.gz -O z -o $output_dir/${analysis_name}_filtNAlt1.vcf.gz
    # Overwrite original vcf with filtered vcf
    mv $output_dir/${analysis_name}_filtNAlt1.vcf.gz $output_dir/${analysis_name}.vcf.gz
    # Remove intermediate file
    rm -f $output_dir/${analysis_name}_filtNAlt1.vcf.gz
fi
## 
# Copy pre-made vcf 
### cp /gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/allpops_r10000_Mono_N12_bootstrap/allpops_r10000_Mono_N12_A99/allpops_r10000_Mono_N12_A99.vcf.gz $output_dir/${analysis_name}.vcf.gz

## Chosen vcf (used to swap out vcfs in bug testing)
vcf_SFS=$output_dir/${analysis_name}

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

# Create projection preview file
# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz --unfolded -p $output_dir/pop_file.txt --preview > $output_dir/proj_$foldtype.txt
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt --preview > $output_dir/proj_$foldtype.txt
fi

# Display projection preview
cat $output_dir/proj_${$foldtype}.txt

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
# runs code but cancels if projection takes longer than 30 seconds. 
# jointMAF are produced quickly but MSFS files can take a long time to produce for large datasets, which are not needed here.
# There is no hope of making a multiSFS with this number of populations, it would be huge (~32 million cells)

# Unfolded or folded
if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
timeout 900s \
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt --unfolded -v -a -f --total-length $SNPcount -o $output_dir/SFS_$foldtype/ --prefix ${analysis_name} \
--proj=$(awk 'NR==1 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==6 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==7 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==8 {print $2}' $output_dir/best_proj_$foldtype.txt)
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
timeout 900s \
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/pop_file.txt -v -a -f --total-length $SNPcount -o $output_dir/SFS_$foldtype/ --prefix ${analysis_name} \
--proj=$(awk 'NR==1 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==2 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==3 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==4 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==5 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==6 {print $2}' $output_dir/best_proj_$foldtype.txt),\
$(awk 'NR==7 {print $2}' $output_dir/best_proj_$foldtype.txt),$(awk 'NR==8 {print $2}' $output_dir/best_proj_$foldtype.txt)
fi


# Deactivate easySFS environment
conda deactivate
module purge

##### 
mkdir -p $output_dir/fsc_run
cd $output_dir/fsc_run

## Remove old obs
rm -f ./*.obs

# Copy over jointMAF file to fsc run directory
cp $output_dir/SFS_$foldtype/fastsimcoal2/${analysis_name}_joint*.obs ./

## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/fsc_run/${analysis_name}.tpl
echo "8 samples to simulate :" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/fsc_run/${analysis_name}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print $1"$"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Samples sizes and samples age" >> $output_dir/fsc_run/${analysis_name}.tpl
awk '{print $2}' $output_dir/best_proj_$foldtype.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/fsc_run/${analysis_name}.tpl
awk '{print "0"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "0" >> $output_dir/fsc_run/${analysis_name}.tpl
# Historical events
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "7 historical event" >> $output_dir/fsc_run/${analysis_name}.tpl
# Get population numbers for events (remember fsc starts counting at 0)
# CLAC merges into LUIB
echo "TDivRWest@ $(awk '$1=="CLAC" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="LUIB" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE1 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
#OBSE merges into DUIN
echo "TDivREast@ $(awk '$1=="OBSE" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE2 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
# CLAM merges into LUIM
echo "TDivMWest@ $(awk '$1=="CLAM" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="LUIM" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE3 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
# OBSM merges into DUIM
echo "TDivMEast@ $(awk '$1=="OBSM" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIM" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE4 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
# Resi west and east merge
echo "TDivResi@ $(awk '$1=="LUIB" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE5 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
# Migration west and east merge
echo "TDivMigr@ $(awk '$1=="LUIM" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIM" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE6 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs@ $(awk '$1=="DUIN" {print NR-1}' $output_dir/pop_uniq.txt) $(awk '$1=="DUIM" {print NR-1}' $output_dir/pop_uniq.txt) 1 RESIZE7 0 0" >> $output_dir/fsc_run/${analysis_name}.tpl

echo "//Number of independent loci [chromosome]" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1 0" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "1" >> $output_dir/fsc_run/${analysis_name}.tpl
echo "//per Block:data typ" >> $output_dir/fsc_run/${analysis_name}.tpl

echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_dir/fsc_run/${analysis_name}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=1000
maxNPOP=1000000

echo "// Priors and rules file" > $output_dir/fsc_run/${analysis_name}.est
echo "// *********************" >> $output_dir/fsc_run/${analysis_name}.est
echo "[PARAMETERS]" >> $output_dir/fsc_run/${analysis_name}.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/fsc_run/${analysis_name}.est
echo "//all N are in number of haploid individuals" >> $output_dir/fsc_run/${analysis_name}.est
## All popuulations
awk -v mxNPOP=$maxNPOP -v mnNPOP=$minNPOP '{print "1 "$1"$ unif "mnNPOP" "mxNPOP" output"}' $output_dir/pop_uniq.txt >> $output_dir/fsc_run/${analysis_name}.est
# East and west ancestral sizes
echo "1 ResiWest$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 ResiEast$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 MigrWest$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 MigrEast$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
# Ecotype ancestral sizes
echo "1 Resi$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 Migr$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 Ancs$ unif $minNPOP $maxNPOP output" >> $output_dir/fsc_run/${analysis_name}.est

# Divergence times (ordered from oldest to most recent)
# The lower range limit is an absolute minimum, whereas the upper range is only used as a
# maximum for choosing a random initial value for this parameter. There is actually no upper
# limit to the search range, as this limit can grow by 30% after each cycle
echo "1 TDivAncs@ unif 100 200000 output" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivMigr@ unif 100 TDivAncs@ output paramInRange" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivResi@ unif 100 TDivAncs@ output paramInRange" >> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivMEast@ unif 100 TDivMigr@ output paramInRange">> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivMWest@ unif 100 TDivMigr@ output paramInRange">> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivREast@ unif 100 TDivResi@ output paramInRange">> $output_dir/fsc_run/${analysis_name}.est
echo "1 TDivRWest@ unif 100 TDivResi@ output paramInRange">> $output_dir/fsc_run/${analysis_name}.est

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE1 = ResiWest$/LUIB$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE2 = ResiEast$/DUIN$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE3 = MigrWest$/LUIM$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE4 = MigrEast$/DUIM$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE5 = Resi$/ResiEast$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE6 = Migr$/MigrEast$ hide" >> $output_dir/fsc_run/${analysis_name}.est
echo "0 RESIZE7 = Ancs$/Migr$ hide" >> $output_dir/fsc_run/${analysis_name}.est


## Run fsc
if [[ $foldtype == "folded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 --foldedSFS -m -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_dir/fsc_run/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ $foldtype == "unfolded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 -d -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_dir/fsc_run/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

############################
 ##### Plot results #####
############################

# Move into fsc run directory
cd $output_dir/fsc_run/${analysis_name}

# Load R module
module load R-uoneasy/4.2.1-foss-2022a
## Plot maxPar file
Rscript ~/code/Github/stickleback-Uist-species-pairs/Helper_scripts/ParFileViewer.R ${analysis_name}_maxL.par $output_dir/pop_uniq.txt

# Copy results to aggregate results directory
mkdir -p $wkdir/results/$vcf_ver/demographic/fastsimcoal2/results_plots/${analysis_dir}/

cp $output_dir/fsc_run/${analysis_name}/${analysis_name}.bestlhoods $wkdir/results/$vcf_ver/demographic/fastsimcoal2/results_plots/${analysis_dir}/
cp $output_dir/fsc_run/${analysis_name}/${analysis_name}_maxL.par $wkdir/results/$vcf_ver/demographic/fastsimcoal2/results_plots/${analysis_dir}/
cp $output_dir/fsc_run/${analysis_name}/${analysis_name}_maxL.par.pdf $wkdir/results/$vcf_ver/demographic/fastsimcoal2/results_plots/${analysis_dir}/

### Then once all jobs are done, run:
## Rscript ~/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/18.15-fastsimcoal-all-bootstrap-plot.R \
##     $wkdir/results/$vcf_ver/demographic/fastsimcoal2/results_plots/${analysis_dir}/ \
##     ${analysis_name}
Rscript ~/code/Github/stickleback-Uist-species-pairs/2_Results/2_1_Population_genomics/18.15-fastsimcoal-all-bootstrap-plot.R \
/gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/results_plots/allpops_r10000_Mono_folded_N12_bootstrap/ \
allpops_r10000_Mono_folded_N12_A