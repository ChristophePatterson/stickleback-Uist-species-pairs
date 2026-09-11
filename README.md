# Patterson et al. (2026) - Divergence in Stickleback Is Linked to Migration, Not Salinity

This repository contains the code used to generate the analyses presented in the manuscript:

Patterson et al. (2026). Divergence in stickleback is linked to migration, not salinity (in preparation).

## Repository Overview

This project originated as an exploratory whole-genome sequencing (WGS) analysis conducted at the start of my postdoctoral research in the MacColl Lab. The repository was later adapted to contain the analysis pipeline required to reproduce the results presented in the manuscript.

To improve clarity and reproducibility, exploratory analyses not directly relevant to the publication have been removed from the current version of the repository. Earlier developmental code remains accessible through the commit: `ade116831f6e0bf18bddfaeac7b2a7b2170c26c4`

As a result, the repository history contains files that were added and subsequently removed during the development of the project.

## Analysis Workflow

Scripts are numbered to indicate their order of execution. Some numbers may be absent because intermediate scripts were removed during repository cleanup; missing numbers do not indicate missing steps in the published workflow.

The primary analysis can be reproduced by running the scripts in numerical order.

## Recommended execution order

`00_data_backup.sh` is optional and is not an analysis step.

### 1. Mapping, depth and variant calling

Run these in order:

1. `01_Sequence_data_download.sh`
2. `02-fastqc.sh`
3. `03-map-reads.sh`
4. `04-read-depth.sh`
5. `04-summary-read-depth.sh` (this invokes `04a-plot-read-depth-summary.R`)
6. `05-clean-bam.sh`
7. `06-read-depth-clean.sh`
8. `06-summary-clean-read-depth.sh` (clean-depth summary and sex-determination plots)
9. `07-snp-calling.sh`
10. `08-concat-filter-snps.sh`
11. `08-vcf2GenomicsGeneral.sh`

The R scripts `08a-plot-individual-variant-stats.R` and `08b-convert-vcf-to-geno.R`
are helper steps invoked by later shell scripts; they are not an additional
numbered stage. `08b-convert-vcf-to-geno.R` is used by `09-SNP-analysis.sh`.

### 2. Population-genomic analyses

Run `09-SNP-analysis.sh` after the filtered VCFs and Genomics General files
exist. The following are parallel branches whose inputs are produced by stages
8–9:

- Sliding-window differentiation: `10-sliding-window-population-distance.sh` -> `10.1-sliding-window-plot.R`
- All-population combinations: `10b-sliding-window-all-population-combinations.sh` -> `10.3-sliding-window-all-pop-combn.R`
- Private alleles: `10c-private-alleles.sh` -> `10.5-Private-alleles-plot.R`
- Sliding-window PCA: `10d-sliding-window-pca.sh` -> `10d1-calculate-sliding-window-pca.R` -> `10d2-plot-sliding-window-pca.R`
- Population heterozygosity: `10e-sliding-window-population-heterozygosity.sh` -> `10e1-plot-population-heterozygosity.R` -> `10e2-plot-adaptive-divergence.R`
- Twisst: `11-twisst-sliding-windows.sh` or `11b-twisst-population-combinations.sh` -> `11a1-summarise-twisst-results.R` -> `11b2-plot-combined-twisst-populations.R`
- NJ sliding-window plots: `11.5-njtree_sliding_window.R` (requires the corresponding window results)

`11-twisst-sliding-windows.sh` and `11b-twisst-population-combinations.sh` are alternative Twisst
workflows, not consecutive mandatory steps.

### 3. Divergence, phylogeny and functional follow-up

After the required 10.x/11.x results are available:

1. CSS: `15.0-CSS.sh` -> `15.1-CSS_plot.R`
2. CSS population sensitivity: `15.2-CSS-Populations.sh` and
   `15.3-CSS-dropPopulations.sh` -> `15.3-CSS-dropPopulations-combine.R` ->
   `15.4-CSS-annotation.R`
3. RAxML-NG: `16.0-RAxML.sh` -> `16.1-RAxML_plot.R`
4. Gene/variant investigation: `17-gene-variant-investigation.sh` ->
   `17a-go-enrichment-analysis.R`
5. Demographic modelling: run `18.0-fastsimcoal-setup.sh`, then the relevant
   SFS-generation branch (`18.0-fastsimcoal-lochs.sh`,
   `18a-fastsimcoal-loch-allsites-sfs.sh`, or
   `18.1.1-fastsimcoal-loch-ecotype-recent.sh`), followed by the selected
   model scripts `18.1`–`18.6` and their plot scripts.
6. Inversion follow-up: `19-chrI-inversion-divergence.sh` ->
   `19.1-ChrI-inv-divergence-est.R`.

`2_Results/2_1_Population_genomics/20-combine-analysis-figures.R` reads products
from several 10.x analyses and should be run as a final figure-combination
step, after those products exist (and, if its panels are included, after the
relevant 15.x–19.x products).

`LiftOver_conversion.sh` and `Helper_scripts/LiftOver-Jones-and-Robert-kingsman.sh`
are reference-coordinate conversion utilities. Run them before any analysis
that consumes their converted coordinates, not as unconditional pipeline
steps.




## Data Availability

Raw sequence files are currenlty being process for upload to Genbank.

## Contact

If you encounter difficulties reproducing the analyses or have questions about specific scripts, please contact:

Christophe Patterson Research Fellow, Population Genomics, University of Nottingham

