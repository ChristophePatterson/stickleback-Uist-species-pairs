#!/bin/bash
# Laura Dean and Christophe Patterson
# 21/13/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5g
#SBATCH --time=01:00:00
#SBATCH --job-name=BD_readdepth_sum
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.2.1-foss-2022a

# set variables
in_filepath=(~/data/sticklebacks/bams/bamstats/QC/raw_bams)
out_filepath=(~/data/sticklebacks/bams/bamstats/QC/raw_bams/Multi-Bam-QC/)
rm -r $out_filepath
mkdir -p $out_filepath

# Define the bigdata file
bigdata="/gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/bigdata_Christophe_2025-03-28.csv"

## Create input config for qualimap
## Uses find to locate all reports and them 
find $in_filepath -wholename *qualimapReport.html | awk -F '/' -v filepath="$in_filepath" '{ print $11 " " filepath "/" $11 "/" }' > qualimap.tmp.txt

## Create summary file
echo -e "sample\treads\tmapped_reads\tpercentage_mapped\tmn_coverage\tstd_coverage\tAve_map_qc\tdupl_reads" > $out_filepath/global_raw_report_custom.txt
cat qualimap.tmp.txt | while read line 
do
    ## Extract filename and filepath
    filepath=$(echo $line | awk '{ print $2 }')
    filename=$(echo $line | awk '{ print $1 }')
    ## Get bam stats file
    QC=$(echo ${filepath}genome_results.txt)
    ## Extract varibles  of interest (any new added make sure to add colname at top)
    # number of reads
    reads=$(sed -n 's/     number of reads = //p' $QC)
    # Number of mapped reads
    map_reads=$(sed -n 's/     number of mapped reads = //p'  $QC)
    # Mean coverage
    mn_cov=$(sed -n 's/     mean coverageData = //p'  $QC)
    ## Standard deviation coverage
    std_cov=$(sed -n 's/     std coverageData = //p'  $QC)
    # Mean mapping quality
    map_qlty=$(sed -n 's/     mean mapping quality = //p'  $QC)
    # Number of duplicated reads
    dup_reads=$(sed -n 's/     number of duplicated reads (estimated) = //p'  $QC)
    ## Print out to one file
    echo -e "${filename}\t${reads}\t${map_reads}\t${mn_cov}\t${std_cov}\t${map_qlty}\t${dup_reads}" >> $out_filepath/global_raw_report_custom.txt
done

Rscript /gpfs01/home/mbzcp2/code/Github/stickleback-Uist-species-pairs/04.2-read-depth-summary-plots.R $out_filepath/global_raw_report_custom.txt $out_filepath


## Run Qualimap
~/apps/qualimap_v2.3/qualimap multi-bamqc -d qualimap.tmp.txt \
    -outformat HTML -outdir $out_filepath -outfile global_raw_report.html


