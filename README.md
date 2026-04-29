# TCGA and METABRIC hypoxia, survival, and CNA analysis

This repository contains the code and derived analysis files used to reproduce the main analyses for the manuscript:

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