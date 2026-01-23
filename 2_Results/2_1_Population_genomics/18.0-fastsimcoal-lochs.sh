#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=60g
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#SBATCH --job-name=fastsimcoal2-ind-lochs.sh
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
model_name=lochs_mono_${foldtype}_nCDS_nHFst_r${randSNP}

## Output
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/$model_name)
mkdir -p $output_dir

# Output directory for all model results
model_output_dir=$output_dir/models_${model_name}_all_results
mkdir -p $model_output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

# Due to overlap in writing files, if slurm array is not equal to 1 then wait 30 seconds
if [ ! $SLURM_ARRAY_TASK_ID = "1" ]; then
   sleep 30
fi

## Get list of populations and samples
## Get unique combination of populations
if [ $SLURM_ARRAY_TASK_ID = "1" ]; then
    ## Get each unique Population (but replacing st with fw)
    grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
        awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $10, $13 } $13=="st" { print $1, $9, "fw" }' > $output_dir/pop_file.txt
    # Get unique waterbodies
    awk '{ print $2 }' $output_dir/pop_file.txt | sort | uniq > $output_dir/waterbody_uniq.txt
    ## Uniq pop names of each ecotype
    # For residential
    awk '$4=="resi" { print $3 }' $output_dir/pop_file.txt | sort | uniq > $output_dir/resi_uniq.txt
    # For anadromous/migratory
    awk '$4=="anad" { print $3 }' $output_dir/pop_file.txt | sort | uniq > $output_dir/anad_uniq.txt
fi

## Get population that equals slurm array
pop=$(awk -v slurmA=$SLURM_ARRAY_TASK_ID 'NR==slurmA {print $0}' $output_dir/waterbody_uniq.txt)

mkdir -p $output_dir/$pop
## Subset to samples
grep -w -f $wkdir/vcfs/$vcf_ver/${species}_samples.txt $output_dir/pop_file.txt | 
   awk -F '\t' -v OFS='\t' -v pop=$pop '{ print $1, $2, $3, $4 }' | \
   sed 's/TOST/TORM/g' | sed 's/OLST/OLAV/g' | awk -v pop=$pop 'NR!=1 && $2==pop {print $1, $3}' | \
   sort -k 2 > $output_dir/$pop/${pop}_pop_file.txt

# Create file with list of individuals
awk '{print $1}' $output_dir/$pop/${pop}_pop_file.txt > $output_dir/$pop/${pop}_ind_file.txt

## Get unique population as variable
pop0=$(awk -v slurmA=$SLURM_ARRAY_TASK_ID 'NR==slurmA { print $1 }' $output_dir/resi_uniq.txt)
pop1=$(awk -v slurmA=$SLURM_ARRAY_TASK_ID 'NR==slurmA { print $1 }' $output_dir/anad_uniq.txt)

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites that don't have a at least some (non-zero) minor allele freq, must filter to just snps first.

## Filter out windows that are in the top 5% of Fst
Fst_calcs=$wkdir/results/$vcf_ver/sliding-window/indPops/sliding_window_w25kb_s5kb_m1_PopPair_APARX.csv
## Figure out which column is the Fst for the two populations
Fst_col=$(head -1 $Fst_calcs | tr ',' '\n' | nl -w1 -s: | grep -w "Fst_${pop0}_${pop1}" | cut -d: -f1)
# Caculate 95th percentile of Fst values
Fst_95th=$(awk -F ',' -v OFS='\t' 'NR>1 { print $'$Fst_col' }' $Fst_calcs | sort -n | awk ' { a[i++]=$1; } END { print a[int(0.95*i)]; }' )

# Subset file to scaffold, start, end, and pop pair Fst and filter out those above 95th percentile
awk -F ',' -v OFS='\t' 'NR!=1 { print $1, $2, $3, $'$Fst_col' }' $Fst_calcs | 
    awk -v Fst_95th=$Fst_95th '$4>=Fst_95th' > $output_dir/$pop/${pop}_fst_filtered_windows.txt

# Remove sites that are
    # 1) Have more than 2 alternate alleles
    # 2) Are in the top 5% of Fst between the two populations
    # 3) Are in coding sequences
    # 4) Randomly prune to set number of bases within a window SNPs
bcftools view -i 'N_ALT<=1' -S $output_dir/$pop/${pop}_ind_file.txt $vcf | \
    bcftools filter -T ^$output_dir/$pop/${pop}_fst_filtered_windows.txt | \
    bcftools filter -T ^/gpfs01/home/mbzcp2/data/sticklebacks/genomes/Duke_GAcu_1_CDS.bed | \
    bcftools +fill-tags -- -t AN,AC,AF,MAF | \
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_dir/$pop/${pop}_r${randSNP}.vcf.gz

# If foldtype is unfolded filter out all alt fixed sites
if [[ ${foldtype} == "unfolded" ]]; then
    bcftools view -e 'N_ALT=1 && AC=AN' $output_dir/${analysis_name}.vcf.gz -O z -o $output_dir/${analysis_name}_filtNAlt1.vcf.gz
    # Overwrite original vcf with filtered vcf
    mv $output_dir/${analysis_name}_filtNAlt1.vcf.gz $output_dir/${analysis_name}.vcf.gz
    # Remove intermediate file
    rm -f $output_dir/${analysis_name}_filtNAlt1.vcf.gz
fi

### cp /gpfs01/home/mbzcp2/data/sticklebacks/results/GCA_046562415.1_Duke_GAcu_1.0_genomic/ploidy_aware_HWEPops_MQ10_BQ20/demographic/fastsimcoal2/lochs_mono_${foldtype}_nCDS_nHFst_r10000/$pop/${pop}_r${randSNP}.vcf.gz \
### $output_dir/$pop/${pop}_r${randSNP}.vcf.gz


## Chosen vcf (used to swap out vcfs in bug testing)
vcf_SFS=$output_dir/$pop/${pop}_r${randSNP}

SAMPcount=$(bcftools query -l $vcf_SFS.vcf.gz | wc -l | awk '{print $1}')
SEQcount=$(bcftools query -l $vcf_SFS.vcf.gz | wc -l | awk '{print $1*2}')

# Calculate number of SNPs input to SFS
SNPcount=$(bcftools view -H $vcf_SFS.vcf.gz | wc -l)

# Deactivate bcftools enviroment
conda deactivate
module purge

conda activate easySFS-env

######################################
##### Get best proj for easySFS #####
######################################

if [[ ${foldtype} == "unfolded" ]]; then
echo "Fold type is: ${foldtype}"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj_unfolded.txt
fi

if [[ ${foldtype} == "folded" ]]; then
echo "Fold type is: ${foldtype}"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj_folded.txt
fi

cat $output_dir/$pop/${pop}_proj_${foldtype}.txt

## Convert projected SFS to long format table 1 for population 0 and 1
grep -A 1 $pop0 $output_dir/$pop/${pop}_proj_${foldtype}.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop0}_proj_long_${foldtype}.txt
grep -A 1 $pop1 $output_dir/$pop/${pop}_proj_${foldtype}.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop1}_proj_long_${foldtype}.txt

# Extract best projection number
topprojpop0=$(awk '{print $2}' $output_dir/$pop/${pop0}_proj_long_${foldtype}.txt | sort -n | tail -1)
topprojpop1=$(awk '{print $2}' $output_dir/$pop/${pop1}_proj_long_${foldtype}.txt | sort -n | tail -1)
## Get best proj for each population
bestprojpop0=$(awk -v topproj=$topprojpop0 '$2==topproj {print $1 }' $output_dir/$pop/${pop0}_proj_long_${foldtype}.txt | sort -n | tail -1)
bestprojpop1=$(awk -v topproj=$topprojpop1 '$2==topproj {print $1 }' $output_dir/$pop/${pop1}_proj_long_${foldtype}.txt | sort -n | tail -1)

SEQpop0=$(tail -n 1 $output_dir/$pop/${pop0}_proj_long_${foldtype}.txt | awk '{print $1}')
SEQpop1=$(tail -n 1 $output_dir/$pop/${pop1}_proj_long_${foldtype}.txt | awk '{print $1}')

echo "Best projection for $pop0 is $bestprojpop0"
echo "Best projection for $pop1 is $bestprojpop1"

############################
##### Run easySFS #####
############################

## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt -a -f --total-length $SNPcount -o $output_dir/$pop/SFS_${foldtype}/ --prefix ${pop0}-${pop1}-$foldtype --proj=$bestprojpop0,$bestprojpop1
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt -a -f --total-length $SNPcount -o $output_dir/$pop/SFS_${foldtype}/ --prefix ${pop0}-${pop1}-$foldtype --proj=$bestprojpop0,$bestprojpop1
fi

# Deactivate easySFS conda
conda deactivate

##########################################################
  ## Create model files for three difference models ##
##########################################################

## Model one -- Isolation (Iso)
## Model two -- Isolation with Migration (IsoMig)
## Model three -- Secondary Contact (IsoSC)

##########################################################
    ### Model 1: Isolation (Iso) ###
##########################################################

## Create Isolation with Migration model files
mkdir -p $output_dir/$pop/Iso
cd $output_dir/$pop/Iso

## Copy MSFS obs file and rename to match model

## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS, needs to be renamed to DSFS for fsc
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-Iso-${foldtype}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-Iso-${foldtype}_MSFS.obs
fi

## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "2 samples to simulate :" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "NPOP0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "NPOP1" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Samples sizes and samples age" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "$bestprojpop0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "$bestprojpop1" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "1" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Migration matrix 0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "2 historical event" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "TDIV 0 1 1 RESIZE 0 0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "TDIV 0 0 0 0 0 0 // kills pop0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Number of independent loci [chromosome]" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "1 0" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "1" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "//per Block:data typ" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl
echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.tpl

## Create estimates file
echo "// Priors and rules file" > $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "// *********************" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "[PARAMETERS]" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "//all N are in number of haploid individuals" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "1 ANCSIZE unif 1000 1000000 output" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "1 NPOP0 unif 1000 1000000 output" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "1 NPOP1 unif 1000 1000000 output" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "1 TDIV unif 1000 2000000 output" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "[COMPLEX PARAMETERS]" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est
echo "0 RESIZE = ANCSIZE/NPOP1 hide" >> $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}.est

echo "Running model: Isolation (Iso) for populations $pop0 and $pop1 with ${foldtype} SFS"
## Run fsc
if [[ ${foldtype} == "folded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-Iso-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-Iso-${foldtype}.est -y 4 --foldedSFS -u -m -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

if [[ ${foldtype} == "unfolded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-Iso-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-Iso-${foldtype}.est -y 4 -u -d -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

## Copy results into single directory
cp $output_dir/$pop/Iso/${pop0}-${pop1}-Iso-${foldtype}/* $model_output_dir


##########################################################
    ### Model 2: Isolation with Migration (IsoMig) ###
##########################################################

## Create Isolation with Migration model files
mkdir -p $output_dir/$pop/IsoMig
cd $output_dir/$pop/IsoMig

## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-IsoMig-${foldtype}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-IsoMig-${foldtype}_MSFS.obs
fi

## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "2 samples to simulate :" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "NPOP0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "NPOP1" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Samples sizes and samples age" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "$bestprojpop0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "$bestprojpop1" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "2" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Migration matrix 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Migration matrix 1" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "0 MIG10" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "MIG01 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "3 historical event" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "TMIG 0 0 1 1 0 1" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "TDIV 0 1 1 RESIZE 0 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "TDIV 0 0 0 0 0 0 // kills pop0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Number of independent loci [chromosome]" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "1 0" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "1" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "//per Block:data typ" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl
echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.tpl

## Create estimates file
echo "// Priors and rules file" > $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "// *********************" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "[PARAMETERS]" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "//all N are in number of haploid individuals" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 ANCSIZE unif 1000 1000000 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 NPOP0 unif 1000 1000000 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 NPOP1 unif 1000 1000000 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 N0M10 logunif 1 NPOP0 hide" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 N1M01 logunif 1 NPOP1 hide" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 TDIV unif 1000 2000000 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "0 TPRP unif 0.5 0.9 hide" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "[COMPLEX PARAMETERS]" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "0 RESIZE = ANCSIZE/NPOP1 hide" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "0 MIG10 = N0M10/NPOP0 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "0 MIG01 = N1M01/NPOP1 output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est
echo "1 TMIG = TDIV*TPRP output" >> $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}.est

echo "Running model: Isolation with Migration (IsoMig) for populations $pop0 and $pop1 with ${foldtype} SFS"
## Run fsc
if [[ ${foldtype} == "folded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-IsoMig-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-IsoMig-${foldtype}.est -y 4 --foldedSFS -u -m -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ ${foldtype} == "unfolded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-IsoMig-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-IsoMig-${foldtype}.est -y 4 -u -d -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

## Copy results into single directory
cp $output_dir/$pop/IsoMig/${pop0}-${pop1}-IsoMig-${foldtype}/* $model_output_dir

##########################################################
    ### Model 3: Secondary Contact (IsoSC) ###
##########################################################

## Create Isolation with Migration model files
mkdir -p $output_dir/$pop/IsoSC
cd $output_dir/$pop/IsoSC

## Create unfolded output
if [[ ${foldtype} == "unfolded" ]]; then
# Create SFS
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-IsoSC-${foldtype}_DSFS.obs
fi

## Create folded output
if [[ ${foldtype} == "folded" ]]; then
# Create SFS
cp $output_dir/$pop/SFS_${foldtype}/fastsimcoal2/${pop0}-${pop1}-${foldtype}_MSFS.obs ./${pop0}-${pop1}-IsoSC-${foldtype}_MSFS.obs
fi

## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "2 samples to simulate :" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "NPOP0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "NPOP1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Samples sizes and samples age" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "$bestprojpop0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "$bestprojpop1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "2" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Migration matrix 0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "0 MIG10" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "MIG01 0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Migration matrix 1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "0 0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "3 historical event" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "TMIG 0 0 1 1 0 1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "TDIV 0 1 1 RESIZE 0 1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "TDIV 0 0 0 0 0 1 // kills pop0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Number of independent loci [chromosome]" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "1 0" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "1" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "//per Block:data typ" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl
echo "FREQ 1 0 5.11e-9 OUTEXP"  >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.tpl

## Create estimates file
echo "// Priors and rules file" > $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "// *********************" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "[PARAMETERS]" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "//all N are in number of haploid individuals" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 ANCSIZE unif 1000 1000000 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 NPOP0 unif 1000 1000000 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 NPOP1 unif 1000 1000000 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 N0M10 logunif 1 NPOP0 hide" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 N1M01 logunif 1 NPOP1 hide" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 TDIV unif 1000 2000000 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "0 TPRP unif 0.1 0.5 hide" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "[COMPLEX PARAMETERS]" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "0 RESIZE = ANCSIZE/NPOP1 hide" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "0 MIG10 = N0M10/NPOP0 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "0 MIG01 = N1M01/NPOP1 output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est
echo "1 TMIG = TDIV*TPRP output" >> $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}.est


echo "Running model: Isolation with Secondary Contact (IsoSC) for populations $pop0 and $pop1 with ${foldtype} SFS"

## Run fsc
if [[ ${foldtype} == "folded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-IsoSC-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-IsoSC-${foldtype}.est -y 4 --foldedSFS -u -m -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi
if [[ ${foldtype} == "unfolded" ]]; then
echo "Fold type is: ${foldtype}"
~/apps/fsc28_linux64/fsc28 -t ${pop0}-${pop1}-IsoSC-${foldtype}.tpl -n 100000 -e ${pop0}-${pop1}-IsoSC-${foldtype}.est -y 4 -u -d -M -L 50 -c $SLURM_CPUS_PER_TASK > fsc_${pop0}-${pop1}-$foldtype_log_jobID${SLURM_ARRAY_TASK_ID}.txt
fi

## Copy results into single directory
cp $output_dir/$pop/IsoSC/${pop0}-${pop1}-IsoSC-${foldtype}/* $model_output_dir