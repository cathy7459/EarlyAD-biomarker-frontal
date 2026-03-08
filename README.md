## Pipeline overview

![Pipeline](https://raw.githubusercontent.com/cathy7459/EarlyAD-biomarker-frontal/main/figures/pipeline.jpg)
# Early-stage Alzheimer’s biomarker discovery using integrative transcriptomic modeling

This repository provides a **reproducible bioinformatics pipeline** for identifying robust early-stage Alzheimer’s disease (AD) biomarkers using integrative analysis of multiple transcriptomic datasets.

The workflow integrates:

- Differential expression analysis
- PCA-guided gene expansion
- Network filtering using WGCNA
- Elastic Net feature selection
- Cross-platform validation using single-nucleus RNA sequencing

The goal is to identify **cross-platform consensus genes consistently detected across independent transcriptomic datasets.**

---

# Overview

Early detection of Alzheimer’s disease remains a major challenge.

This project aims to identify **robust transcriptomic biomarkers associated with early-stage AD progression** by integrating multiple GEO datasets and validating findings across different platforms.

The pipeline consists of three main stages:

1. Bulk transcriptomic biomarker discovery  
2. Machine-learning based feature selection  
3. Cross-platform validation using single-nucleus RNA-seq  

---

# Pipeline

## 1. Bulk transcriptomic biomarker discovery

Dataset: **GSE118553**

### Differential expression analysis

Three contrasts were tested using **limma**:

- AsymAD vs Control  
- AD vs Control  
- AD vs AsymAD  

### Preclinical Core Gene Set (PCGS)

The preclinical core genes were defined as:

```
PCGS = DEG(AsymAD vs Control)
       ∩ (DEG(AD vs Control) ∪ DEG(AD vs AsymAD))
```

These genes represent **early transcriptional changes associated with AD progression.**

### PCA-based dimensionality reduction

Principal component analysis was applied to identify major transcriptional axes.

The PCA step serves two purposes:

1. Dimensionality reduction  
2. Extraction of **PC loadings representing major biological variation**

### PCA-guided gene expansion

Genes highly correlated with selected PCs were added to the feature space.

This step expands the PCGS to capture genes **co-regulated with early AD transcriptional programs.**

### Network filtering using WGCNA

The expanded gene set was filtered using **Weighted Gene Co-expression Network Analysis (WGCNA)**.

This step identifies **co-expression modules associated with disease progression** and removes noise genes.

The resulting set is referred to as the **Core Gene Network**.

---

## 2. Independent dataset projection and machine learning

Training datasets:

- GSE33000  
- GSE132903  

Steps:

1. Projection of candidate genes onto independent bulk datasets  
2. Probe-to-gene mapping and gene-level expression collapse  
3. Data standardization using training statistics  
4. Feature selection using Elastic Net logistic regression  
5. Model training and performance evaluation  

Elastic Net was used instead of LASSO to better account for correlated gene expression features.

### External validation

The trained model was evaluated on:

- **GSE122063**

Model performance was measured using **ROC curves and AUC scores**.

---

## 3. Cross-platform validation using single-nucleus RNA-seq

Datasets:

- GSE188545  
- GSE243292  

Steps:

1. Projection of bulk-derived hub genes onto snRNA-seq datasets  
2. Hub gene filtering to minimize technical noise  
3. **LOOCV-based hub score evaluation**  
4. Independent external validation using **GSE5281**

---

## 4. Final biomarker identification

Genes detected consistently across bulk transcriptomic datasets and single-cell validation were defined as:

**Cross-platform consensus genes**

These genes represent **robust early-stage Alzheimer’s disease biomarkers.**

---

# Repository structure

```
scripts/
    01_pcgs/             preclinical core gene set
    02_pca_expansion/    PCA-guided gene expansion
    03_wgcna/            co-expression network filtering
    04_elastic_net/      Elastic Net modeling
    05_cross_platform/   snRNA cross-platform validation

R/utils/                 helper functions
config/                  configuration parameters
data/                    raw and processed datasets
outputs/                 tables and figures
```

---

# Requirements

R ≥ 4.2

Required packages:

- GEOquery  
- limma  
- WGCNA  
- glmnet  
- pROC  
- data.table  
- matrixStats  
- dplyr  
- tibble  
- stringr  

---

# Installation example

```r
install.packages(c(
  "data.table",
  "dplyr",
  "tibble",
  "stringr",
  "matrixStats"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "WGCNA"
))
```

---

# Running the pipeline

```
Rscript scripts/run_all.R
```

Outputs will be written to:

```
outputs/
```

---

# Key findings

- The **temporal cortex** showed stronger diagnostic signals than the frontal cortex.  
- **3 frontal hub genes** and **9 temporal hub genes** were consistently detected across independent datasets.  
- External validation achieved **AUC up to 0.94**.  
- Cross-platform validation confirmed robust biomarkers across bulk and single-cell transcriptomic datasets.

---

# Reproducibility

All parameters are defined in:

```
config/config.R
```

Raw data are not stored in this repository.

GEO datasets can be downloaded from:

https://www.ncbi.nlm.nih.gov/geo/
