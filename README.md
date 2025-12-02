# Code for: Multi-Modal Profiling of Alzheimer’s Disease: From Dual Molecular Signatures to Neuropathology-Predominant Cognitive Prediction 

This repository contains the complete R and Python code used to perform the analyses and generate the figures for the manuscript "Multi-Modal Profiling of Alzheimer’s Disease From Dual Molecular Signatures to Neuropathology-Predominant Cognitive Prediction".

## Overview

This study presents an integrative multi-modal framework for characterizing and prioritizing biomarkers in Alzheimer's Disease (AD). Our analysis combines single-nucleus transcriptomics, quantitative neuropathology, spatial transcriptomics, and clinical metadata. The code in this repository covers four main analytical modules:
1.  Cellular composition overview
2.  Pseudobulk differential expression and functional enrichment analysis
3.  Donor-level machine learning for predictive modeling
4.  Statistical analysis of spatial transcriptomic data

## Data Availability

The raw snRNA-seq and spatial transcriptomics data are publicly available from the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) and associated public repositories, as detailed in the manuscript's Methods section.

**All processed data files required to run the analysis scripts in this repository are available on Zenodo:**
* **DOI:** 10.5281/zenodo.17299986

## Software and Package Versions

All analyses were conducted using the following software and package versions to ensure reproducibility.

### R Environment
* **R version 4.5.1** (2025-06-13)
* **Key Packages:**
    * brms `v2.23.0`
    * DESeq2 `v1.48.2`
    * dplyr `v1.1.4`
    * ggplot2 `v4.0.0`
    * Seurat `v5.3.0`
    * spdep `v1.4.1`
    * *Additional packages listed in the scripts include: `Matrix v1.7.4`, `HDF5Array v1.36.0`, `SummarizedExperiment v1.38.1`, `tibble v3.3.0`, `tidyr v1.3.1`, `stringr v1.5.2`, `purrr v1.1.0`, `lme4 v1.1.37`, `lmerTest v3.1.3`, `rstan v2.32.7`, `sf v1.0.21`, `spatstat v3.4.0`, `patchwork v1.3.2`.*

### Python Environment
* **Python version 3.9.6**
* **Key Packages:**
    * numpy `v1.26.2`
    * pandas `v2.3.2`
    * scanpy `v1.10.3`
    * scikit-learn (`sklearn`) `v1.6.1`
    * scipy `v1.13.1`
    * matplotlib `v3.8.2`

## Instructions for Use / Workflow

The scripts are organized by analysis module and should be run in the following general order. Please ensure the processed data from the Zenodo repository has been downloaded and placed in the appropriate directory as referenced in the scripts.

**1. Cell Composition Analysis:**
* Run `Omics proportion(Bar Chart).R` and `Omics proportion(Pie Chart).R` to generate the overview of cellular composition (Figure 1).

**2. Differential Expression Analysis:**
* Run `DEG&GO preprocessing.ipynb` to pre-processing the required data for R or you can also directly use the procossed data offered on Zenodo.
* Run `DESeq2 Differential Analysis.R` to perform the differential expression and functional enrichment analysis (Figure 2, Figure 3, Supplementary Tables S1-S5).

**3. Machine Learning Analysis:**
* Run `Machine Learning （GroupKFold).ipynb` to perform the donor-level machine learning benchmark and incremental value analysis.

**4. Spatial Analysis:**
* Run `Spatial Transcriptomics Preprocessing.ipynb` to pre-processing the required data for R or you can also directly use the procossed data offered on Zenodo.
* Run `Spatial related processing.R` to perform the Beta GLMM for layer occupancy and the Moran's I analysis (Figure 4, Table 2).

## Citation

If you use this code or the associated data in your research, please cite our paper:

> "Author(s), (Year). Title of Paper. *Journal Name*. DOI"
