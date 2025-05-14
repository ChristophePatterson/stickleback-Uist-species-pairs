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
#SBATCH --job-name=TTcalculator
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out


tabix /gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/Obse_347/Obse_347_chr1.vcf.gz
tabix /gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/Obsm_641/Obsm_641_chr1.vcf.gz

mkdir -p $wkdir/results/TTmethod/vcfs/test

bcftools merge \
-o $wkdir/results/TTmethod/vcfs/test/OBSE_OBSM_chr1.vcf.gz \
-O v \
--threads $SLURM_CPUS_PER_TASK \
/gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/Obse_347/Obse_347_chr1.vcf.gz \
/gpfs01/home/mbzcp2/data/sticklebacks/results/TTmethod/vcfs/Obsm_641/Obsm_641_chr1.vcf.gz

tabix $wkdir/results/TTmethod/vcfs/test/OBSE_OBSM_chr1.vcf.gz

mkdir -p $wkdir/results/TTmethod/vcfs/test/sfs
mkdir -p $wkdir/results/TTmethod/vcfs/test/TTcals
start_pos=(1)
inc_pos=(100000)
end_pos=$((start_pos+inc_pos))
# stop_pos=$(bcftools query  -f '%POS\n' $wkdir/results/TTmethod/vcfs/test/OBSE_OBSM_chr1.vcf.gz | tail -n 1)
echo $end_pos
echo $stop_pos

while [ "$end_pos" -le "$stop_pos" -o "$start_pos" -ge "$stop_pos" ]
do
    bcftools view -s Obse_347,Obsm_641 -r NC_053212.1:$start_pos-$end_pos $wkdir/results/TTmethod/vcfs/test/OBSE_OBSM_chr1.vcf.gz | bcftools view --max-alleles 2 --exclude-types indels -e 'GT=="./."'| \
    bcftools query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > $wkdir/results/TTmethod/vcfs/test/sfs/OBSE_OBSM_chr1_${start_pos}-${end_pos}.sfs

    # Rscriot tt.R <path to unfolded 2dsfs> <prefix for output files> <assumed mutation rate> <assumed generation time>
    Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/2_Results/2_2_Divergence_estimates/TT-calculator.R \
    $wkdir/results/TTmethod/vcfs/test/sfs/OBSE_OBSM_chr1_${start_pos}-${end_pos}.sfs \
    $wkdir/results/TTmethod/vcfs/test/TTcals/OBSE_OBSM_chr1_${start_pos}-${end_pos} \
    "2.1e-08" \
    "1"

    end_pos=$((end_pos+inc_pos))
    start_pos=$((start_pos+inc_pos))
    echo $start_pos

done


