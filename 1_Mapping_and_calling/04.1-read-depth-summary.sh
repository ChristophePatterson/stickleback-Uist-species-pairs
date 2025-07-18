#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_readdepth_sum
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.2.1-foss-2022a

# set variables
# genome_name=(GCF_016920845.1_GAculeatus_UGA_version5_genomic)
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)

in_filepath=(~/data/sticklebacks/bams/$genome_name/bamstats/QC/raw_bams)
out_filepath=(~/data/sticklebacks/bams/$genome_name/bamstats/QC/raw_bams/Multi-Bam-QC)
rm -r $out_filepath
mkdir -p $out_filepath

## Create input config for qualimap
# Output header
echo -e "sample\tbam_file\treads\tmapped_reads\tpercentage_mapped\tmn_coverage\tstd_coverage\tAve_map_qc\tdupl_reads" > "$out_filepath/global_raw_report_custom.txt"

pairdata=(~/data/sticklebacks/bams/$genome_name/species_pairs_sequence_data.csv)

# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

while read line; do
# Extract individual using awk
    individual=$(echo $line | awk -F ',' '{ print $1 }')
    echo "Processing: $individual"
    QC="${in_filepath}/$individual/genome_results.txt"
    if [[ ! -f "$QC" ]]; then
        echo "QC file not found: $QC"
        continue
    fi
    echo $QC
    bam_file=$(grep -m1 "bam file =" "$QC" | awk -F' = ' '{print $2}')
    reads=$(grep -m1 "number of reads =" "$QC" | awk -F' = ' '{print $2}')
    map_reads=$(grep -m1 "number of mapped reads =" "$QC" | awk -F' = ' '{print $2}')
    mn_cov=$(grep -m1 "mean coverageData =" "$QC" | awk -F' = ' '{print $2}')
    std_cov=$(grep -m1 "std coverageData =" "$QC" | awk -F' = ' '{print $2}')
    map_qlty=$(grep -m1 "mean mapping quality =" "$QC" | awk -F' = ' '{print $2}')
    dup_reads=$(grep -m1 "number of duplicated reads (estimated) =" "$QC" | awk -F' = ' '{print $2}')

    echo -e "${individual}\t${bam_file}\t${reads}\t${map_reads}\t\t${mn_cov}\t${std_cov}\t${map_qlty}\t${dup_reads}" >> "$out_filepath/global_raw_report_custom.txt"
done < $pairdata

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/1_Mapping_and_calling/04.2-read-depth-summary-plots.R $out_filepath/global_raw_report_custom.txt $out_filepath

awk -v OFS='\t' -v inpath="$in_filepath" ' NR!=1 { print $1, inpath"/"$1"/" }' $out_filepath/global_raw_report_custom.txt > $out_filepath/qualimap.tmp.txt
## Run Qualimap
~/apps/qualimap_v2.3/qualimap multi-bamqc -d $out_filepath/qualimap.tmp.txt --feature-file $genome.bed \
    -outformat HTML -outdir $out_filepath -outfile global_raw_report.html


