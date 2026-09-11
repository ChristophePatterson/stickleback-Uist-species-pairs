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

`00.00-backup-data.sh` is optional and is not an analysis step.

### 1. Mapping, depth and variant calling

Run these in order:

1. `01.00-download-sequence-data.sh`
2. `02.00-run-fastqc.sh`
3. `03.00-map-reads.sh`
4. `04.00-calculate-read-depth.sh`
5. `04.01-summarise-read-depth.sh` (this invokes `04.02-plot-read-depth-summary.R`)
6. `05.00-clean-bam.sh`
7. `06.00-calculate-clean-read-depth.sh`
8. `06.01-summarise-clean-read-depth.sh` (clean-depth summary and sex-determination plots)
9. `07.00-call-snps.sh`
10. `08.00-concatenate-and-filter-snps.sh`
11. `08.03-convert-vcf-to-genomics-general.sh`

The R scripts `08.01-plot-individual-variant-stats.R` and `08.02-convert-vcf-to-geno.R`
are helper steps invoked by later shell scripts; they are not an additional
numbered stage. `08.02-convert-vcf-to-geno.R` is used by `09.00-run-snp-analysis.sh`.

### 2. Population-genomic analyses

Run `09.00-run-snp-analysis.sh` after the filtered VCFs and Genomics General files
exist. The following are parallel branches whose inputs are produced by stages
8–9:

- Sliding-window differentiation: `10.00-sliding-window-population-distance.sh` -> `10.01-plot-sliding-window-population-distance.R`
- All-population combinations: `10.10-sliding-window-population-combinations.sh` -> `10.11-plot-sliding-window-population-combinations.R`
- Private alleles: `10.20-private-alleles.sh` -> `10.21-plot-private-alleles.R`
- Sliding-window PCA: `10.30-sliding-window-pca.sh` -> `10.31-calculate-sliding-window-pca.R` -> `10.32-plot-sliding-window-pca.R`
- Population heterozygosity: `10.40-sliding-window-population-heterozygosity.sh` -> `10.41-plot-population-heterozygosity.R` -> `10.42-plot-adaptive-divergence.R`
- Twisst: `11.00-twisst-sliding-windows.sh` or `11.10-twisst-population-combinations.sh` -> `11.01-summarise-twisst-results.R` -> `11.11-plot-combined-twisst-populations.R`
- NJ sliding-window plots: `11.20-plot-njtree-sliding-windows.R` (requires the corresponding window results)

`11.00-twisst-sliding-windows.sh` and `11.10-twisst-population-combinations.sh` are alternative Twisst
workflows, not consecutive mandatory steps.

### 3. Divergence, phylogeny and functional follow-up

After the required `10.xx`/`11.xx` results are available:

1. CSS: `15.00-calculate-css.sh` -> `15.01-plot-css.R`
2. CSS population sensitivity: `15.10-calculate-css-populations.sh` and
   `15.20-calculate-css-drop-populations.sh` -> `15.21-combine-css-drop-populations.R` ->
   `15.22-annotate-css-regions.R`
3. RAxML-NG: `16.00-run-raxml.sh` -> `16.01-plot-raxml-results.R`
4. Gene/variant investigation: `17.00-investigate-gene-variants.sh` ->
   `17.01-run-go-enrichment-analysis.R`
5. Demographic modelling: run `18.00-run-loch-pair-model-selection.sh`
   with an SFS type and loch name, once for each combination of
   `folded`/`unfolded` and `CLAC`/`DUIN`/`LUIB`/`OBSE`. For example:
   `sbatch 18.00-run-loch-pair-model-selection.sh folded CLAC`.
   This script compares the five demographic models using only the resident
   and migratory ecotypes within each loch.
6. Plot model selection with `18.01-plot-loch-pair-model-selection.R`,
   passing the corresponding `models_loch_Mselect_*_all_results` directory
   as its only argument.
7. Inversion follow-up: `19.00-extract-chrI-inversion-variants.sh` ->
   `19.01-analyse-chrI-inversion-divergence.R`.

`2_Results/2_1_Population_genomics/20.00-combine-analysis-figures.R` reads products
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
