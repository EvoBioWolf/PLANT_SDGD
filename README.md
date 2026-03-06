This repository contains the analysis pipeline for the paper "Species richness buffers losses of genetic diversity in experimental grassland communities".

The workflow processes raw Illumina sequencing data (mcSCRB-seq), performs variant calling using a customized pipeline, and analyzes species-genetic diversity correlations (SGDC) using data from the Jena Experiment.

## Repository Structure & Workflow

The scripts are numbered `01` through `08` to follow the chronological order of the analysis pipeline described in the *Methods* section of the manuscript.

### 1. Pre-processing & Demultiplexing

- `01_demultiplex.sh`: Splits raw sequencing data into sample-specific files based on barcodes.

### 2. Read Processing & Deduplication

*Corresponds to "Read processing and variant calling" in Methods.*

- `02_extractUMI_50bp.py` / `.sh`: Because libraries target transcript 3' ends, this step extracts the Unique Molecular Identifiers (UMIs) and handles the asymmetric paired-end reads (where Read 1 contains the barcode/UMI and Read 2 contains the transcript).
- `02_cdhit_dedup.sh`: Removes PCR duplicates using `cd-hit-dup`, utilizing the UMIs extracted in the previous step.
- `02_getR2.sh`: Retains only Read 2 (the transcript sequence) for mapping, removing potential poly(A) stretches.

### 3. Alignment & Expression Quantification

- `03_mappingToRef.sh`: Maps cleaned reads to species-specific PacBio Iso-Seq reference transcriptomes using BWA.
- `03_countMapping.sh` / `countOverlapped.py`: Quantifies transcript expression levels (used later as a covariate in GLMMs).

### 4. Variant Calling

*Corresponds to GATK workflow described in Methods.*

- `04_individualCalling.sh`: Runs GATK HaplotypeCaller in GVCF mode for each individual sample.
- `04_cohortCalling.sh`: Combines GVCFs using GATK GenotypeGVCFs to create a multi-sample VCF for each species.

### 5. Quality Control

- `05_replicatesTest.[py/R/sh]`: Assesses genotype call repeatability using technical replicates. Species with insufficient repeatability were excluded from further analysis here.

### 6. Heterozygosity Analysis (GLMM)

*Corresponds to "Statistical analyses" in Methods.*

- `06_withinSp_heterozygosity.py`: Calculates the number of heterozygous sites per transcript/individual from the VCF files.
- `06_withinSp_heterozygosity_fitModel.R`: Fits the Zero-Inflated Negative Binomial GLMMs (`glmmTMB`) to test the effect of species richness, biomass, and expression on heterozygosity.
- `06_withinSp_heterozygosity_metaFor.R`: Performs the random-effects meta-analysis to estimate the pooled mean effect of species richness across all species.
- `06_withinSp_heterozygosity_Fig2a.R`: Generates Figure 2A (Slope estimates).

### 7. Genetic Differentiation ($F_{ST}$)

*Corresponds to "Genetic differentiation analysis" in Methods.*

- `07_plotLevel_Fst.py`: Calculates pairwise $F_{ST}$ between experimental plots based on allele frequencies.
- `07_plotLevel_Fst_Figs.R`: Generates Figure 2B ($F_{ST}$ comparisons between low and high diversity plots).

### 8. Plot-Level Genetic Diversity (Expected Heterozygosity) 

*Corresponds to "Plot-level expected heterozygosity" in Methods.* 

- `08_plotLevel_expected_heterozygosity.py`: Calculates the mean expected heterozygosity ($H_{exp}$) for each experimental plot.     

- `08_plotLevel_expected_heterozygosity.sh`: Bash runner to execute the Python script across multiple species VCFs. 

- `08_plotLevel_expected_heterozygosity_submit.sh`: HPC cluster submission script.



### 9. Simulation (Drift Modeling)

*Corresponds to "Ne simulation analysis" in Methods.*

- `09_simulation.R`: Simulates the expected loss of genetic diversity due to the "dilution effect" (reduced population size) under neutral drift, and compares empirical data to these theoretical curves (Figure 3).

## Metadata Files

- `plot_div.txt`: Contains metadata for the experimental plots (e.g., Species Richness level (1, 2, 4, 8, 16, 60), Plot IDs).
- `sp_functionalGroup.txt`: Classification of the 19 focal species (Grasses, Legumes, Herbs).
- `Figure2_species_order_factor.rds`: R data object preserving the ordering of species for consistent plotting.

## Dependencies

The following software versions were used in this analysis:

- BWA (0.7.17-r1188)

- GATK (4.3.0)

- CD-HIT (4.8.1)

- Cutadapt (4.6)

- R

  (4.2.2)

  - Key packages: `glmmTMB`, `metafor`, `tidyverse`

- Python (3.x)

## Usage Note on `_submit.sh` Files

Files ending in `_submit.sh` (e.g., `04_cohortCalling_submit.sh`) are wrapper scripts used to submit jobs to a High-Performance Computing (HPC) cluster (e.g., via Slurm or PBS). The core logic is contained in the corresponding `.sh`, `.R`, or `.py` files.
