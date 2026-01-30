#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10g
#SBATCH --time=24:00:00
#SBATCH --array=1-100
#SBATCH --job-name=fastsimcoal2-CLDO-CDLO-COLD-sigMig
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
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/model_selection_MSFS_${foldtype}_r${randSNP})

output_model=($output_dir/model_files/run_${foldtype}_A${SLURM_ARRAY_TASK_ID})

## Check is SFS directory exists
if [ ! -d $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/ ]; then
    echo "SFS directory does not exist: $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/"
    echo "Please run 18.0-fastsimcoal-setup.sh first to generate the SFS."
    exit 1
fi
## Create models with all three combinations of relationships between the four resident populations

## Models
# 1) mCLDO
# 2) mCDLO
# 3) mCOLD

#############################
    ### Model 1) mCLDO ###
#############################

analysis_name=SFS_${SLURM_ARRAY_TASK_ID}_sigMig-mCLDO_${foldtype}

mkdir -p $output_model
cd $output_model

## Copy MSFS obs file and rename to match model
## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS, needs to be renamed to DSFS for fsc
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_MSFS.obs
fi


echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_model/${analysis_name}.tpl
echo "5 samples to simulate :" >> $output_model/${analysis_name}.tpl
echo "//Population effective sizes (number of genes)" >> $output_model/${analysis_name}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print $1"$"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Samples sizes and samples age" >> $output_model/${analysis_name}.tpl
awk '{print $2}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_best_proj_long_unfolded.txt >> $output_model/${analysis_name}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_model/${analysis_name}.tpl
awk '{print "0"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_model/${analysis_name}.tpl
echo "0" >> $output_model/${analysis_name}.tpl
# Historical events
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_model/${analysis_name}.tpl
echo "4 historical event" >> $output_model/${analysis_name}.tpl
# Get population numbers for events (remember fsc starts counting at 0)
# CLAC merges into LUIB
echo "TDivCL@ $(awk '$1=="CLAC" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="LUIB" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE1 0 0" >> $output_model/${analysis_name}.tpl
#OBSE merges into DUIN
echo "TDivDO@ $(awk '$1=="OBSE" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE2 0 0" >> $output_model/${analysis_name}.tpl
# Resi west and east merge
echo "TDivResi@ $(awk '$1=="LUIB" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE3 0 0" >> $output_model/${analysis_name}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs@ $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="anad" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE4 0 0" >> $output_model/${analysis_name}.tpl

echo "//Number of independent loci [chromosome]" >> $output_model/${analysis_name}.tpl
echo "1 0" >> $output_model/${analysis_name}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_model/${analysis_name}.tpl
echo "1" >> $output_model/${analysis_name}.tpl
echo "//per Block:data typ" >> $output_model/${analysis_name}.tpl

echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_model/${analysis_name}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=100
maxNPOP=10000000

echo "// Priors and rules file" > $output_model/${analysis_name}.est
echo "// *********************" >> $output_model/${analysis_name}.est
echo "[PARAMETERS]" >> $output_model/${analysis_name}.est
echo "//#isInt? #name #dist.#min #max" >> $output_model/${analysis_name}.est
echo "//all N are in number of haploid individuals" >> $output_model/${analysis_name}.est
# Ecotype ancestral sizes
echo "1 Ancs$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 Resi$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
# East and west ancestral sizes
echo "1 ResiCL$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 ResiDO$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est

## All popuulations
awk -v mxNPOP=$maxNPOP '{print "1 "$1"$ unif 1000 "mxNPOP" output"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.est

# Divergence times
# The lower range limit is an absolute minimum, whereas the upper range is only used as a
# maximum for choosing a random initial value for this parameter. There is actually no upper
# limit to the search range, as this limit can grow by 30% after each cycle
echo "1 TDivAncs@ unif 1000 200000 output" >> $output_model/${analysis_name}.est
echo "1 TDivResi@ unif 1000 TDivAncs@ output paramInRange" >> $output_model/${analysis_name}.est
echo "1 TDivCL@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est
echo "1 TDivDO@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_model/${analysis_name}.est
echo "0 RESIZE1 = ResiCL$/LUIB$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE2 = ResiDO$/DUIN$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE3 = Resi$/ResiDO$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE4 = Ancs$/anad$ hide" >> $output_model/${analysis_name}.est

################################
 ####  Run fsc  ####
################################

if [[ $foldtype == "folded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 --foldedSFS --multiSFS -m -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ $foldtype == "unfolded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 -d --multiSFS -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

############################
 ##### Plot results #####
############################

# Move into fsc run directory
cd $output_model/${analysis_name}

# Load R module
module load R-uoneasy/4.2.1-foss-2022a
## Plot maxPar file
Rscript ~/code/Github/stickleback-Uist-species-pairs/Helper_scripts/ParFileViewer.R ${analysis_name}_maxL.par $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt

# Copy results to aggregate results directory
mkdir -p $output_dir/results_plots/

cp ${analysis_name}.bestlhoods $output_dir/results_plots/
cp ${analysis_name}_maxL.par $output_dir/results_plots/
cp ${analysis_name}_maxL.par.pdf $output_dir/results_plots/


#############################
    ### Model 2) mCLDO ###
#############################

analysis_name=SFS_${SLURM_ARRAY_TASK_ID}_sigMig-mCDLO_${foldtype}

mkdir -p $output_model
cd $output_model

## Copy MSFS obs file and rename to match model
## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS, needs to be renamed to DSFS for fsc
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_MSFS.obs
fi

echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_model/${analysis_name}.tpl
echo "5 samples to simulate :" >> $output_model/${analysis_name}.tpl
echo "//Population effective sizes (number of genes)" >> $output_model/${analysis_name}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print $1"$"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Samples sizes and samples age" >> $output_model/${analysis_name}.tpl
awk '{print $2}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_best_proj_long_unfolded.txt >> $output_model/${analysis_name}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_model/${analysis_name}.tpl
awk '{print "0"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_model/${analysis_name}.tpl
echo "0" >> $output_model/${analysis_name}.tpl
# Historical events
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_model/${analysis_name}.tpl
echo "4 historical event" >> $output_model/${analysis_name}.tpl
# Get population numbers for events (remember fsc starts counting at 0)
# CLAC merges into LUIB
echo "TDivCD@ $(awk '$1=="CLAC" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE1 0 0" >> $output_model/${analysis_name}.tpl
#OBSE merges into DUIN
echo "TDivLO@ $(awk '$1=="LUIB" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="OBSE" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE2 0 0" >> $output_model/${analysis_name}.tpl
# Resi west and east merge
echo "TDivResi@ $(awk '$1=="OBSE" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE3 0 0" >> $output_model/${analysis_name}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs@ $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="anad" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE4 0 0" >> $output_model/${analysis_name}.tpl

echo "//Number of independent loci [chromosome]" >> $output_model/${analysis_name}.tpl
echo "1 0" >> $output_model/${analysis_name}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_model/${analysis_name}.tpl
echo "1" >> $output_model/${analysis_name}.tpl
echo "//per Block:data typ" >> $output_model/${analysis_name}.tpl

echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_model/${analysis_name}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=100
maxNPOP=10000000

echo "// Priors and rules file" > $output_model/${analysis_name}.est
echo "// *********************" >> $output_model/${analysis_name}.est
echo "[PARAMETERS]" >> $output_model/${analysis_name}.est
echo "//#isInt? #name #dist.#min #max" >> $output_model/${analysis_name}.est
echo "//all N are in number of haploid individuals" >> $output_model/${analysis_name}.est
# Ecotype ancestral sizes
echo "1 Ancs$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 Resi$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
# East and west ancestral sizes
echo "1 ResiCD$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 ResiLO$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est

## All popuulations
awk -v mxNPOP=$maxNPOP '{print "1 "$1"$ unif 1000 "mxNPOP" output"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.est

# Divergence times
# The lower range limit is an absolute minimum, whereas the upper range is only used as a
# maximum for choosing a random initial value for this parameter. There is actually no upper
# limit to the search range, as this limit can grow by 30% after each cycle
echo "1 TDivAncs@ unif 1000 200000 output" >> $output_model/${analysis_name}.est
echo "1 TDivResi@ unif 1000 TDivAncs@ output paramInRange" >> $output_model/${analysis_name}.est
echo "1 TDivCD@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est
echo "1 TDivLO@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_model/${analysis_name}.est
echo "0 RESIZE1 = ResiCD$/DUIN$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE2 = ResiLO$/OBSE$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE3 = Resi$/ResiLO$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE4 = Ancs$/anad$ hide" >> $output_model/${analysis_name}.est

################################
 ####  Run fsc  ####
################################

if [[ $foldtype == "folded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 --foldedSFS --multiSFS -m -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ $foldtype == "unfolded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 -d --multiSFS -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

############################
 ##### Plot results #####
############################

# Move into fsc run directory
cd $output_model/${analysis_name}

# Load R module
module load R-uoneasy/4.2.1-foss-2022a
## Plot maxPar file
Rscript ~/code/Github/stickleback-Uist-species-pairs/Helper_scripts/ParFileViewer.R ${analysis_name}_maxL.par $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt

# Copy results to aggregate results directory
mkdir -p $output_dir/results_plots/

cp ${analysis_name}.bestlhoods $output_dir/results_plots/
cp ${analysis_name}_maxL.par $output_dir/results_plots/
cp ${analysis_name}_maxL.par.pdf $output_dir/results_plots/

#############################
    ### Model 3) mCOLD ###
#############################

analysis_name=SFS_${SLURM_ARRAY_TASK_ID}_sigMig-mCOLD_${foldtype}

mkdir -p $output_model
cd $output_model

## Copy MSFS obs file and rename to match model
## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS, needs to be renamed to DSFS for fsc
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_sigMig_$foldtype/fastsimcoal2/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_${foldtype}_MSFS.obs $output_model/${analysis_name}_MSFS.obs
fi

echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_model/${analysis_name}.tpl
echo "5 samples to simulate :" >> $output_model/${analysis_name}.tpl
echo "//Population effective sizes (number of genes)" >> $output_model/${analysis_name}.tpl
## Print out NPOP lines from waterbody uniq file
awk '{print $1"$"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Samples sizes and samples age" >> $output_model/${analysis_name}.tpl
awk '{print $2}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/SFS_${SLURM_ARRAY_TASK_ID}_sigMig_best_proj_long_unfolded.txt >> $output_model/${analysis_name}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_model/${analysis_name}.tpl
awk '{print "0"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_model/${analysis_name}.tpl
echo "0" >> $output_model/${analysis_name}.tpl
# Historical events
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_model/${analysis_name}.tpl
echo "4 historical event" >> $output_model/${analysis_name}.tpl
# Get population numbers for events (remember fsc starts counting at 0)
# CLAC merges into LUIB
echo "TDivCO@ $(awk '$1=="CLAC" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="OBSE" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE1 0 0" >> $output_model/${analysis_name}.tpl
#OBSE merges into DUIN
echo "TDivLD@ $(awk '$1=="LUIB" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE2 0 0" >> $output_model/${analysis_name}.tpl
# Resi west and east merge
echo "TDivResi@ $(awk '$1=="OBSE" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE3 0 0" >> $output_model/${analysis_name}.tpl
# Resi and Migr merge into Ancestral
echo "TDivAncs@ $(awk '$1=="DUIN" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) $(awk '$1=="anad" {print NR-1}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt) 1 RESIZE4 0 0" >> $output_model/${analysis_name}.tpl

echo "//Number of independent loci [chromosome]" >> $output_model/${analysis_name}.tpl
echo "1 0" >> $output_model/${analysis_name}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_model/${analysis_name}.tpl
echo "1" >> $output_model/${analysis_name}.tpl
echo "//per Block:data typ" >> $output_model/${analysis_name}.tpl

echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_model/${analysis_name}.tpl

############################
 ## Create estimates file ##
############################

# Set min and max population sizes
minNPOP=100
maxNPOP=10000000

echo "// Priors and rules file" > $output_model/${analysis_name}.est
echo "// *********************" >> $output_model/${analysis_name}.est
echo "[PARAMETERS]" >> $output_model/${analysis_name}.est
echo "//#isInt? #name #dist.#min #max" >> $output_model/${analysis_name}.est
echo "//all N are in number of haploid individuals" >> $output_model/${analysis_name}.est
# Ecotype ancestral sizes
echo "1 Ancs$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 Resi$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
# East and west ancestral sizes
echo "1 ResiCO$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est
echo "1 ResiLD$ unif $minNPOP $maxNPOP output" >> $output_model/${analysis_name}.est

## All popuulations
awk -v mxNPOP=$maxNPOP '{print "1 "$1"$ unif 1000 "mxNPOP" output"}' $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt >> $output_model/${analysis_name}.est

# Divergence times
# The lower range limit is an absolute minimum, whereas the upper range is only used as a
# maximum for choosing a random initial value for this parameter. There is actually no upper
# limit to the search range, as this limit can grow by 30% after each cycle
echo "1 TDivAncs@ unif 1000 200000 output" >> $output_model/${analysis_name}.est
echo "1 TDivResi@ unif 1000 TDivAncs@ output paramInRange" >> $output_model/${analysis_name}.est
echo "1 TDivCO@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est
echo "1 TDivLD@ unif 1000 TDivResi@ output paramInRange">> $output_model/${analysis_name}.est

#  Complex parameters
echo "[COMPLEX PARAMETERS]" >> $output_model/${analysis_name}.est
echo "0 RESIZE1 = ResiCO$/OBSE$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE2 = ResiLD$/DUIN$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE3 = Resi$/ResiLD$ hide" >> $output_model/${analysis_name}.est
echo "0 RESIZE4 = Ancs$/anad$ hide" >> $output_model/${analysis_name}.est

################################
 ####  Run fsc  ####
################################

if [[ $foldtype == "folded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 --foldedSFS -m --multiSFS -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ $foldtype == "unfolded" ]]; then
~/apps/fsc28_linux64/fsc28 -t ${analysis_name}.tpl -n 100000 -e ${analysis_name}.est -y 4 -d --multiSFS -M -L 50 -c $SLURM_CPUS_PER_TASK > $output_model/fsc_${analysis_name}_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

############################
 ##### Plot results #####
############################

# Move into fsc run directory
cd $output_model/${analysis_name}

# Load R module
module load R-uoneasy/4.2.1-foss-2022a
## Plot maxPar file
Rscript ~/code/Github/stickleback-Uist-species-pairs/Helper_scripts/ParFileViewer.R ${analysis_name}_maxL.par $output_dir/SFS/SFS_${SLURM_ARRAY_TASK_ID}/pop_sigMig_uniq.txt

# Copy results to aggregate results directory
mkdir -p $output_dir/results_plots/

cp ${analysis_name}.bestlhoods $output_dir/results_plots/
cp ${analysis_name}_maxL.par $output_dir/results_plots/
cp ${analysis_name}_maxL.par.pdf $output_dir/results_plots/