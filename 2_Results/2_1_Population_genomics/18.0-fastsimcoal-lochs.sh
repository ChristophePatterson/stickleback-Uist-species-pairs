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
#SBATCH --job-name=fastsimcoal2-setup
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

## Output
output_dir=($wkdir/results/$vcf_ver/demographic/fastsimcoal2/lochs_mono_r${randSNP})
mkdir -p $output_dir

## Input vcf
vcf=$wkdir/vcfs/$vcf_ver/stickleback_all.NOGTDP5.MEANGTDP5_200.Q60.SAMP0.8.MAF2.vcf.gz

## Get list of populations and samples
## Get unique combination of populations
if [ ! -f $output_dir/pop_uniq.txt_combn.txt ]; then
    ## Get each unique Population (but replacing st with fw)
    grep -f $wkdir/vcfs/$vcf_ver/${species}_subset_samples.txt /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-04-28.csv | 
        awk -F ',' -v OFS='\t' '$13!="st" { print $1, $9, $10, $13 } $13=="st" { print $1, $9, "fw" }' > $output_dir/pop_file.txt
    # Get unique waterbodies
    awk '{ print $2 }' $output_dir/pop_file.txt | sort | uniq > $output_dir/waterbody_uniq.txt
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
pop0=$(awk '{ print $2 }' $output_dir/$pop/${pop}_pop_file.txt | sort | uniq | awk 'NR==1 { print $0 }')
pop1=$(awk '{ print $2 }' $output_dir/$pop/${pop}_pop_file.txt | sort | uniq | awk 'NR==2 { print $0 }')

## Subset vcf to specific population
conda activate bcftools-env

# Filter to those specific samples
# Removing sites that don't have a at least some (non-zero) minor allele freq, must filter to just snps first.
# With random filtering for reduced input
bcftools view -i 'N_ALT<=1' -S $output_dir/$pop/${pop}_ind_file.txt $vcf | \
    # bcftools +fill-tags -- -t AN,AC,AF,MAF | \
    # bcftools view -q '0.00000001:minor' -Q '0.9999999:minor' |
    bcftools +prune -n 1 -N rand -w ${randSNP}bp -O z -o $output_dir/$pop/${pop}_r${randSNP}.vcf.gz

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

if [[ $foldtype == "unfolded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj_unfolded.txt
fi

if [[ $foldtype == "folded" ]]; then
echo "Fold type is: $foldtype"
python ~/apps/easySFS/easySFS.py -a -i $vcf_SFS.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt --preview > $output_dir/$pop/${pop}_proj_folded.txt
fi

cat $output_dir/$pop/${pop}_proj_folded.txt

## Convert projected SFS to long format table 1 for population 0 and 1
grep -A 1 $pop0 $output_dir/$pop/${pop}_proj_folded.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop0}_proj_long_$foldtype.txt
grep -A 1 $pop1 $output_dir/$pop/${pop}_proj_folded.txt | grep '(' | sed 's/(/\n/g' | sed 's/)//g' | awk -F ',' '{print $1, $2}' > $output_dir/$pop/${pop1}_proj_long_$foldtype.txt

# Extract best projection number
topprojpop0=$(awk '{print $2}' $output_dir/$pop/${pop0}_proj_long_$foldtype.txt | sort -n | tail -1)
topprojpop1=$(awk '{print $2}' $output_dir/$pop/${pop1}_proj_long_$foldtype.txt | sort -n | tail -1)
## Get best proj for each population
bestprojpop0=$(awk -v topproj=$topprojpop0 '$2==topproj {print $1 }' $output_dir/$pop/${pop0}_proj_long_$foldtype.txt | sort -n | tail -1)
bestprojpop1=$(awk -v topproj=$topprojpop1 '$2==topproj {print $1 }' $output_dir/$pop/${pop1}_proj_long_$foldtype.txt | sort -n | tail -1)

SEQpop0=$(tail -n 1 $output_dir/$pop/${pop0}_proj_long_$foldtype.txt | awk '{print $1}')
SEQpop1=$(tail -n 1 $output_dir/$pop/${pop1}_proj_long_$foldtype.txt | awk '{print $1}')

echo "Best projection for $pop0 is $bestprojpop0"
echo "Best projection for $pop1 is $bestprojpop1"

############################
##### Run easySFS #####
############################

## Create unfolded output
if [[ $foldtype == "unfolded" ]]; then
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz --unfolded -p $output_dir/$pop/${pop}_pop_file.txt -a -f --total-length $SNPcount -o $output_dir/$pop/SFS_$foldtype/ --prefix ${pop}_$foldtype --proj=$bestprojpop0,$bestprojpop1
fi

## Create folded output
if [[ $foldtype == "folded" ]]; then
# Create SFS
python ~/apps/easySFS/easySFS.py -i $vcf_SFS.vcf.gz -p $output_dir/$pop/${pop}_pop_file.txt -a -f --total-length $SNPcount -o $output_dir/$pop/SFS_$foldtype/ --prefix ${pop}_$foldtype --proj=$bestprojpop0,$bestprojpop1
fi

mkdir -p $output_dir/$pop/fsc_run
cd $output_dir/$pop/fsc_run

# Sum of all polymorphic sites included in SFS
# Getting third line
# Cutting first SFS value (monomorphic sites)
# Transforming into awk to sum (%.13 increases decimal places)

## Remove old obs
rm -f ./*.obs
## Create new SFS that has zero in monomorphic sites
#### awk 'NR<3 {print $0}' $output_dir/$pop/SFS_$foldtype/fastsimcoal2/${pop}_${foldtype}_MSFS.obs > ${pop0}-${pop1}_MSFS.obs 
#### echo "0 $(awk 'NR==3 {print $0}' $output_dir/$pop/SFS_$foldtype/fastsimcoal2/${pop}_${foldtype}_MSFS.obs | cut -d ' ' -f 2-)" >> ${pop0}-${pop1}_MSFS.obs 

cp $output_dir/$pop/SFS_$foldtype/fastsimcoal2/${pop}_${foldtype}_MSFS.obs ./${pop0}-${pop1}_MSFS.obs 

# Sum up all values in SFS (non including monomorphic sites)
projSFSsites=$(awk 'NR==3 {print $0}' ${pop0}-${pop1}_MSFS.obs | \
    sed 's/ /\n/g' | \
    awk -v OFMT=%.13g '{sum += $1} END {print sum}')

## Create model parameters file
echo "//Parameters for the coalescence simulation program : simcoal.exe" > $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "2 samples to simulate :" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Population effective sizes (number of genes)" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "NPOP1" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "NPOP2" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Samples sizes and samples age" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "$bestprojpop0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "$bestprojpop1" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Growth rates: negative growth implies population expansion" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Number of migration matrices : 0 implies no migration between demes" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "2" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Migration matrix 0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "0 MIG21" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "MIG12 0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Migration matrix 1" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "0 0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "0 0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "1 historical event" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "TDIV 0 1 1 RESIZE 0 1" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Number of independent loci [chromosome]" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "1 0" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "1" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "//per Block:data typ" >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl
echo "FREQ $projSFSsites 0 5.11e-9 OUTEXP"  >> $output_dir/$pop/fsc_run/$pop0-$pop1.tpl

## Create estimates file
echo "// Priors and rules file" > $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "// *********************" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "[PARAMETERS]" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "//#isInt? #name #dist.#min #max" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "//all N are in number of haploid individuals" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 ANCSIZE unif 1000 10000000 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 NPOP1 unif 1000 10000000 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 NPOP2 unif 1000 10000000 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 N1M21 logunif 1e-2 20 hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 N2M12 logunif 1e-2 20 hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "1 TDIV unif 1000 2000000 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "[COMPLEX PARAMETERS]" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 RESIZE = ANCSIZE/NPOP2 hide" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 MIG21 = N1M21/NPOP1 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est
echo "0 MIG12 = N2M12/NPOP2 output" >> $output_dir/$pop/fsc_run/$pop0-$pop1.est

## Run fsc
~/apps/fsc28_linux64/fsc28 -t $pop0-$pop1.tpl -n 100000 -e $pop0-$pop1.est -y 4 -m -u -M -L 100 -c $SLURM_CPUS_PER_TASK -q 