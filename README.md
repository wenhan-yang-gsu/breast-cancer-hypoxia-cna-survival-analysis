# TCGA and METABRIC hypoxia, survival, and CNA analysis

This repository contains the code and derived analysis files used to reproduce the main analyses:

**Tumor hypoxia is associated with global copy-number alteration burden and subtype-dependent overall survival in breast cancer: evidence from TCGA and METABRIC**

## 1. Overview

This analysis evaluates the association between transcriptome-derived tumor hypoxia, genome-wide copy-number alteration (CNA) burden, and overall survival in breast cancer using TCGA and METABRIC cohorts.

The repository includes:

- the final R analysis script,
- processed per-sample analysis datasets,
- the final model result workbook,
- a plain-text summary of the main results,
- dataset manifest information,
- R session information for reproducibility.

## 2. Required source data

This analysis uses two public breast cancer datasets downloaded in cBioPortal format:

1. **Breast Invasive Carcinoma (TCGA, PanCancer Atlas)**  
   Study identifier: `brca_tcga_pan_can_atlas_2018`

2. **Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)**  
   Study identifier: `brca_metabric`

The full raw cBioPortal study folders are not redistributed in this repository. To rerun the full pipeline from raw source files, place the two cBioPortal study folders inside the same project directory as the R script.

Expected directory structure:

```text
project_directory/
├── brca_tcga_pan_can_atlas_2018/
├── brca_metabric/
├── reproducible_code.R
└── README.md
```

## 3. Repository files

This repository includes the following reproducibility files:

- `reproducible_code.R`: R script used to prepare the analysis datasets, run survival and CNA analyses, and generate the final result files.
- `00_final_results_PHsafe.xlsx`: final model result workbook. The `preferred_results` sheet contains the main estimates used for manuscript reporting.
- `00_final_summary_PHsafe.txt`: plain-text summary of the main final model results.
- `TCGA_processed_survival_dataset.csv`: processed TCGA survival analysis dataset.
- `METABRIC_processed_survival_dataset.csv`: processed METABRIC survival analysis dataset.
- `TCGA_LumB_CNA_processed_dataset.csv`: processed TCGA Luminal B CNA analysis dataset.
- `METABRIC_CNA_processed_dataset.csv`: processed METABRIC CNA analysis dataset.
- `dataset_manifest.csv`: dataset and processed-file manifest.
- `sessionInfo.txt`: R version and package information.

The processed CSV files contain the derived per-sample variables used in the manuscript analyses, including survival variables, hypoxia scores or hypoxia groups, clinical covariates, and CNA burden summaries where applicable.

## 4. How to rerun the analysis

After downloading the two raw cBioPortal study folders and placing them in the expected directory structure, run the following command in R or RStudio:

```r
source("reproducible_code.R")
```
