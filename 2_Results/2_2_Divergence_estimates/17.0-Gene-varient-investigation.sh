#!/bin/bash
# Christophe Patterson
# 27/02/26
# For running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --time=23:00:00
#SBATCH --job-name=Gene_variant_investigation
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out


# load bcftools
module purge
source /gpfs01/home/${USER}/.bashrc
module load singularity/3.8.5
conda activate bcftools-env
module load samtools-uoneasy/1.18-GCC-12.3.0

# Set variables
wkdir=/gpfs01/home/mbzcp2/data/sticklebacks
call_pars=(ploidy_aware_HWEPops_MQ10_BQ20)
genome_name=(GCA_046562415.1_Duke_GAcu_1.0_genomic)
vcf_ver=($genome_name/$call_pars)

CSS_dir="$wkdir/results/$vcf_ver/sliding-window/CSS/dropPops"
output_dir="$wkdir/results/$vcf_ver/Regions_of_interest"

# Draft genome to use
genome_dir=$wkdir/genomes
# genome_name=(GCF_016920845.1_GAculeatus_UGA_version5_genomic)
genome=($genome_dir/$genome_name.fna)

# Make output directory
mkdir -p $output_dir

# Copy over the regions that are significant across all dropped populations
# And convert to .BED format
cat ${CSS_dir}/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.txt | \
  awk -v OFS='\t' -F ',' 'NR != 1 {print $3,$4,$5}' > $output_dir/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.BED

# Create bedgraph for IGV
cat ${CSS_dir}/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.txt | \
  awk -v OFS='\t' -F ',' 'NR != 1 {print $3,$4,$5,$6}' > $output_dir/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.bedgraph

# Use the regions generated from the CSS analysis to extract variants from the vcf file for the regions of interest
# Limit SNPs to those will high variabilty but only samples from DUIN
bcftools view -R $output_dir/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.BED \
  -s Uist22617,Uist22631,Uist22628,Uist22616,Uist22627,Uist22629,Uist22618,Uist22619,Uist22620,Uist22635,Uist22625,Uist22632 \
  -v snps $wkdir/vcfs/$vcf_ver/stickleback.bcf |
  bcftools +fill-tags -- -t AN,AC,AF,MAF | \
  bcftools view -i 'N_ALT<=1' |
  bcftools view --min-ac 6:minor -Oz -o ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_igv.vcf.gz
tabix ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_igv.vcf.gz

# Limit SNPs to those which are high variabilty across all samples
bcftools view -R $output_dir/stickleback.dropPops..wnd2500.sld500.mnSNP1.mthbasepair-mds.MAF0.05_CSS_all_sig_top_regions.BED \
  -v snps $wkdir/vcfs/$vcf_ver/stickleback.bcf |
  bcftools +fill-tags -- -t AN,AC,AF,MAF | \
  bcftools view -i 'N_ALT<=1' |
  bcftools view --min-ac 6:minor -Oz -o ${output_dir}/stickleback_ALL_minAC6_all_sig_top_regions_igv.vcf.gz
tabix ${output_dir}/stickleback_ALL_minAC6_all_sig_top_regions_igv.vcf.gz

## Extract the bam files for the just chrI inversion for the two highest coverage DUIN samples
# Resi sample Uist22628
samtools view -b $wkdir/bams/$genome_name/raw_bams/Uist22628_raw.bam \
  -o ${output_dir}/Uist22628_CM102076.1_26500000-27200000_raw.bam CM102076.1:26500000-27200000
samtools index ${output_dir}/Uist22628_CM102076.1_26500000-27200000_raw.bam
# Migrant sample Uist22617
samtools view -b $wkdir/bams/$genome_name/raw_bams/Uist22617_raw.bam \
  -o ${output_dir}/Uist22617_CM102076.1_26500000-27200000_raw.bam CM102076.1:26500000-27200000
samtools index ${output_dir}/Uist22617_CM102076.1_26500000-27200000_raw.bam
# Migrant sample Uist22631
samtools view -b $wkdir/bams/$genome_name/raw_bams/Uist22631_raw.bam \
  -o ${output_dir}/Uist22631_CM102076.1_26500000-27200000_raw.bam CM102076.1:26500000-27200000
samtools index ${output_dir}/Uist22631_CM102076.1_26500000-27200000_raw.bam


# Convert Fst DUIN AND DUIN to bedgraph for IGV
# Get the column number for Fst_DUIN_DUIM
cat $wkdir/results/$vcf_ver/sliding-window/indPops/sliding_window_w25kb_s5kb_m1_PopPair_APARX.csv |
    awk 'NR == 1 {print}' |
    tr ',' '\n' | grep -n "Fst_DUIN_DUIM"

# Extract the chromosome, start, end and Fst_DUIN_DUIM columns from the sliding window results file, excluding the header and any rows where Fst_DUIN_DUIM is "nan", and save it as a bedgraph file for IGV
awk -F ',' 'NR != 1 && $45 != "nan" {print $1, $2, $3, $45}' \
 $wkdir/results/$vcf_ver/sliding-window/indPops/sliding_window_w25kb_s5kb_m1_PopPair_APARX.csv > ${output_dir}/stickleback_DUIN_DUIM_FST.bedgraph

# Get gtf downloaded from NCBI
# For DUKE this needs some heafty formatting
# Get the chromosome names from the gtf file
awk '{print $1}' $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1.gtf | uniq
# Remove "chr" from chromosome names in the gtf file, and remove leading "0" from chromosome names in the gtf file
awk -F '\t' -v OFS='\t' '{ gsub("chr", "", $1); gsub(/^0+/, "", $1); print }' $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1.gtf > $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames.gtf
# Swap the chromosome names in the gtf file with those in the reference genome
# Get sequencing report for the reference genome
awk 'NR != 1 {print $4}' $wkdir/genomes/GCA_046562415.1/GCA_046562415.1_seq_report.tsv # Print chromosome names from the sequence report to check
awk 'NR != 1 {print $7}' $wkdir/genomes/GCA_046562415.1/GCA_046562415.1_seq_report.tsv # Print chromosome names from the sequence report to check
awk -F '\t' -v OFS='\t' 'NR==FNR{a[$4]=$7; next} {if($1 in a) $1=a[$1]}1' $wkdir/genomes/GCA_046562415.1/GCA_046562415.1_seq_report.tsv $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames.gtf  > $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf
# Check that the chromosome names in the gtf file now match those in the sequence report
awk '{print $1}' $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf | uniq

# Format gtf into zipped file for VEP
grep -v "#" $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf | awk -F '\t' -v OFS='\t' '$3 != "start_codon" && $3 != "intron" {print}' | \
   sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf.gz
tabix $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf.gz

# If its a transcript line, add "transcript_id " to column 9 and put the oringal value in quotation marks
# And remove the trail ".tX" from the transcript to add in the gene_id for VEP
# e.g turn transcript_id 'g3823.t1' into 'transcript_id "g3823.t1"; gene_id "g3823"'
zcat $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed.gtf.gz  |awk -F'\t' -v OFS='\t' '
    $3 == "gene" {
        $9 = "gene_id \"" $9 "\""
    }
    $3 == "transcript" {
        # Add transcript_id
        tid = $9
        $9 = "transcript_id \"" tid "\""

        # Extract gene_id by removing .tX
        gid = gensub(/\.t[0-9]+$/, "", "g", tid)
        $9 = $9 "; gene_id \"" gid "\""
    }
    { print }
' > $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf

bgzip -c $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf > $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf.gz
tabix $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf.gz

# Run VEP
singularity exec -B $output_dir:/output ~/apps/vep.sif vep -i ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_igv.vcf.gz \
  -gtf $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf.gz --fasta $genome --output_file ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions.vep.out --force_overwrite

# Run VEP
singularity exec -B $output_dir:/output ~/apps/vep.sif vep -i ${output_dir}/stickleback_ALL_minAC6_all_sig_top_regions_igv.vcf.gz \
  -gtf $wkdir/genomes/GCA_046562415.1/Duke_GAcu_1_ChrNames_fixed_gene_id.gtf.gz --fasta $genome --output_file ${output_dir}/stickleback_ALL_minAC6_all_sig_top_regions.vep.out --force_overwrite


## Extract missense_variants from the VEP output file, but include headers in the output file
grep -E "^#|missense_variant|stop_gained|stop_lost" ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions.vep.out > ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_coding_consequences.vep.out

# Convert to bed file for IGV
# Extract the chromosome, start, end and gene name from the VEP output file, excluding the header, and save it as a bed file for IGV
# Take second column and split into chromomosome and position, and add 1 to the position to get the end position, and extract the gene name from the INFO column
grep -v "^#" ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_coding_consequences.vep.out | awk -F '\t' -v OFS='\t' '{split($2, a, ":"); print a[1], a[2]-1, a[2], $11}' > ${output_dir}/stickleback_DUIN_minAC6_all_sig_top_regions_coding_consequences.bed
