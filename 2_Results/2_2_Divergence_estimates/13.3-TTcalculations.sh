#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=48:00:00
#SBATCH --array=1,18,31,40,45
#SBATCH --job-name=TTcalculator
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## Load modules
module load R-uoneasy/4.2.1-foss-2022a
module load bcftools-uoneasy/1.18-GCC-13.2.0

## Set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks 
species=stickleback

## Pop Combination file
PopComb=($wkdir/results/TTmethod/vcfs/top_coverage_samples_PopComb.txt)
TopSamps=($wkdir/results/TTmethod/vcfs/top_coverage_samples.txt)
# Reference chromosome
genome_info=(/gpfs01/home/mbzcp2/data/sticklebacks/genomes/GCF_016920845.1_sequence_report.tsv)


## Make array equal to the number of population pairs to make comparison
Pop1=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $PopComb | awk '{ print $1 }')
Pop2=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $PopComb | awk '{ print $2 }')

sample1=$(grep $Pop1 $TopSamps | awk '{ print $1 }')
sample2=$(grep $Pop2 $TopSamps | awk '{ print $1 }')

# sample1=$(awk)
echo "Calculating TT for ${Pop1} and ${Pop2} using $sample1 and $sample2"

# Number of bases to shift up by 
inc_pos=2500000 # Equals 5Mb
# Convert to text in KB
inc_pos_txt=$(expr $inc_pos / 1000)
echo "Running TT on windows of ${inc_pos_txt}Kb"

## Create folder for vcf of each population
vcf_out=($wkdir/results/TTmethod/vcfs/PopComb/${Pop1}_${Pop2})
top_dir=($wkdir/results/TTmethod/TTresults_JK_1wd_${inc_pos_txt}kb)
result_out=($top_dir/${Pop1}_${Pop2})
mkdir -p $vcf_out
mkdir -p $result_out


### Loop throug indivudal chromosomes (exculusing X (19) and Y (22))
for chr in {1..18} {20..21}
do
    # Extract scaffold name from genome info file (need to add on to number fo avoid header)
    scaf=$(awk -F ',' -v line=$((chr + 1)) "FNR==line" $genome_info | awk '{ print $9 }')
    echo "Starting Chromosome $chr using scafffold $scaf"
    # Create text chrI for easier to need file extraction code
    chrI=$(echo chr$chr)

    # If vcf hasn't been indiexed index for sample 1
    if [ ! -f $wkdir/results/TTmethod/vcfs/${sample1}/${sample1}_${chrI}.vcf.gz.tbi ]; then
        	tabix $wkdir/results/TTmethod/vcfs/${sample1}/${sample1}_${chrI}.vcf.gz
    fi
    # and for sample 2
    if [ ! -f $wkdir/results/TTmethod/vcfs/${sample2}/${sample2}_${chrI}.vcf.gz.tbi ]; then
        	tabix $wkdir/results/TTmethod/vcfs/${sample2}/${sample2}_${chrI}.vcf.gz
    fi

    ## Create merged bcf for each of the two samples
    bcftools merge \
        -o $vcf_out/${Pop1}_${Pop2}_${chrI}.vcf.gz \
        -O v \
        --threads $SLURM_CPUS_PER_TASK \
        /gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/${sample1}/${sample1}_${chrI}.vcf.gz \
        /gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/${sample2}/${sample2}_${chrI}.vcf.gz

    tabix $vcf_out/${Pop1}_${Pop2}_${chrI}.vcf.gz

    ## Remove indels and sites that have more than 2 allele varients and sites with missing base calles
    bcftools view --max-alleles 2 --exclude-types indels -e 'GT=="./."' -O v -o $vcf_out/${Pop1}_${Pop2}_filt_${chrI}.vcf.gz $vcf_out/${Pop1}_${Pop2}_${chrI}.vcf.gz
    tabix $vcf_out/${Pop1}_${Pop2}_filt_${chrI}.vcf.gz
    # remove old vcf and tbi
    rm $vcf_out/${Pop1}_${Pop2}_${chrI}.vcf.gz
    rm $vcf_out/${Pop1}_${Pop2}_${chrI}.vcf.gz.tbi

    # Create output directories
    mkdir -p $result_out/sfs
    mkdir -p $result_out/TTcals
    mkdir -p $result_out/TTresults

    # Start position for chromosome
    start_pos=1
    # Starting end position (start position + increment)
    end_pos=$((start_pos+inc_pos))
    # Hight called bases in vcf to stop at
    stop_pos=$(bcftools query  -f '%POS\n' $vcf_out/${Pop1}_${Pop2}_filt_${chrI}.vcf.gz | tail -n 1)
    
    ## Starting at the first section calculate TT for all sections of the chromosome
    while [ "$start_pos" -le "$stop_pos" ]; do
        ## Remove 
        bcftools view -t ${scaf}:$start_pos-$end_pos $vcf_out/${Pop1}_${Pop2}_filt_${chrI}.vcf.gz | \
        bcftools query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > $result_out/sfs/${Pop1}_${Pop2}_${chrI}_${start_pos}-${end_pos}.sfs

        # Rscript TT-calculator.R <path to unfolded 2dsfs> <prefix for output files> <assumed mutation rate> <assumed generation time>
        Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/13.4-TT-calculator.R \
        $result_out/sfs/${Pop1}_${Pop2}_${chrI}_${start_pos}-${end_pos}.sfs \
        $result_out/TTcals/${Pop1}_${Pop2}_${chrI}_${start_pos}-${end_pos} \
        "2.1e-08" \
        "1"
        ## Increase section of chromsome being used
        end_pos=$((end_pos+inc_pos))
        start_pos=$((start_pos+inc_pos))
        # echo $start_pos

    done
## End loop through chromosomes
done

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/13.5-TT-plot.R \
        $result_out \
        $result_out/TTresults/${Pop1}_${Pop2}

## Copy plots into single directory
mkdir -p $top_dir/plots
cp $result_out/TTresults/${Pop1}_${Pop2}* $top_dir/plots


