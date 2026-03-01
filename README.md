[README.md](https://github.com/user-attachments/files/25644238/README.md)# EarlyAD-biomarker-PCA-WGCNA
Reproducible transcriptomic pipeline for early AD biomarker discovery using PCA expansion and WGCNA.
#PCGS README
[README.md](https://github.com/user-attachments/files/25642571/README.md)
# PCGS (Preclinical Core Gene Set) — GSE118553 (Frontal)

This repo contains a cleaned, reproducible R pipeline to:

1. Download **GSE118553** from GEO
2. Keep **frontal** region samples
3. Run **limma** differential expression for 3 contrasts
4. Define a **preclinical core** set:

   `core = DEG(AsymAD vs Control) ∩ (DEG(AD vs Control) ∪ DEG(AD vs AsymAD))`

5. Map **probe → gene symbol** using the platform (GPL) annotation and collapse probes to genes
6. Export core gene expression matrices
7. Expand the feature space using **PCA expansion** (max |cor| with selected PCs)

## Files

- `scripts/run_pcgs_frontal.R` — main pipeline script (edit parameters at the top)
- `R/utils.R` — helper functions
- `outputs/` — default output folder (generated)

## Requirements

R packages:

- GEOquery
- limma
- matrixStats
- dplyr
- stringr
- tibble

Install (example):

```r
install.packages(c("dplyr","stringr","tibble","matrixStats"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma"))
```

## Run

From repo root:

```bash
Rscript scripts/run_pcgs_frontal.R
```

Outputs are written to `outputs/`.

## Key outputs

- `deg_AsymAD_vs_Control.tsv`
- `deg_AD_vs_Control.tsv`
- `deg_AD_vs_AsymAD.tsv`
- `core_probes.txt`, `core_genes.txt`
- `core_expression_gene_x_sample.csv`
- `expanded_genes.txt`
- `expanded_expression_gene_x_sample.csv`
- `pcgs_frontal_pipeline.rds`

## Notes

- Region matching uses a regex (`frontal|pfc|prefrontal|dlpfc|dorsolateral`). If your `pData` uses different naming, change `frontal_regex`.
- The pipeline keeps logic explicit and avoids overwriting objects in-place (one of the biggest readability issues in the original script).

# PCA and WGCNA
#[Uploadin# GSE118553 PCA expansion + WGCNA (cleaned)

This repository is a cleaned, GitHub-friendly refactor of your original script:
`original/PCA_expansion_and_WGCNA_gse118553.R`.

## What this repo does
1. **PCA expansion**: starting from a *core gene set*, it expands to a capped gene set (default 800) using:
   - simple correlations to selected PCs
   - partial correlations (PCk controlling other PCs)
   - AND filter across ≥ `m` PCs, then hard-cap by a combined score

2. **WGCNA**: runs a signed WGCNA on the (expanded) gene set and saves:
   - gene → module mapping
   - module eigengenes
   - a bundled RDS of the full results

## What you must provide upstream
This repo intentionally **does not** re-download GEO data (so it's fast and modular).
Before running, create these objects in your environment:

- `expr_mapped` : probe-level expression matrix (probes x samples)
- `map_df`      : mapping table with column `symbol` (gene symbol) aligned to `expr_mapped` rows
- `expr_gene`   : **core gene** expression matrix (genes x samples), used to fix sample order
- `core_features` : character vector of core gene symbols
- (optional) `group_f` : labels for each sample (names must be sample IDs)

## Run
After your preprocessing script has created the objects above:

```r
source("scripts/run_all.R")
```

Outputs are written to `outputs/`.

## Outputs
- `outputs/expanded_genes.txt`
- `outputs/expanded_expression_gene_x_sample.csv`
- `outputs/pca_expansion_diagnostics.rds`
- `outputs/WGCNA_GeneModuleMap.csv`
- `outputs/WGCNA_ModuleEigengenes.csv`
- `outputs/WGCNA_results.rds`
- (optional) `outputs/WGCNA_ME_trait_correlation.tsv` (+ pvalues)

## Tuning
Edit parameters at the top of `scripts/01_pca_expansion.R` and `scripts/02_wgcna.R`:
- PCA expansion: `K_fixed`, `drop_PC1`, `m_and`, `target_max`
- WGCNA: `minModuleSize`, `mergeCutHeight`, `networkType`

---
If you want, I can also:
- add `optparse` CLI (`--target_max 800` etc.)
- add `renv` to lock package versions for full reproducibility
g README.md…]()

#LASSO
[README.md](https://github.com/user-attachments/files/25659777/README.md)
# LASSO pipeline (GSE33000) — paper-ready scaffold

This repository contains a **reproducible** (paper-style) pipeline to:

1) build a **gene × sample** expression matrix for **GSE33000** from the GEO series matrix  
2) map probes → gene symbols (GPL4372) and **collapse multiple probes per gene**  
3) run **effect-size filtering** (limma; |logFC| + FDR) to define a stable feature set  
4) train a **binomial LASSO** model (glmnet) and (optionally) validate on another GEO dataset

> Raw data are not stored in this repository. Place downloaded files under `data/raw/`.

## Quick start

### 0) Requirements
- R >= 4.2
- Packages: `data.table`, `GEOquery`, `limma`, `glmnet`, `pROC`

### 1) Put inputs in place
Create folders locally (they are gitignored):
- `data/raw/`  
- `data/processed/`  
- `outputs/`

Place:
- `data/raw/GSE33000_series_matrix.txt.gz`
- (optional) a validation series matrix: e.g. `data/raw/GSE44770_series_matrix.txt.gz`
- a core gene list CSV: `data/processed/core_genes.csv` with a column **GeneSymbol**

If you **cannot** reliably infer case/control labels from GEO metadata, create:
- `data/processed/sample_metadata_gse33000.csv` with columns:
  - `sample_id` (e.g., GSMxxxxxxx matching series matrix columns)
  - `condition` (values: `AD` or `Control`)

### 2) Run everything
From the repo root:
```bash
Rscript scripts/run_all.R
```

## Outputs
- `data/processed/gse33000_expr_gene_x_sample.csv`
- `outputs/tables/deg_gse33000.csv`
- `outputs/tables/deg_filtered_genes.csv`
- `outputs/tables/lasso_selected_genes.csv`
- `outputs/figures/roc_train.png` (+ optional validation ROC)

## Reproducibility notes
- All file paths and thresholds are centralized in `config/config.R`.
- The pipeline is split into small scripts under `scripts/` and reusable functions under `R/`.

## Original script
See `original/LASSO with gse33000_1.R` for the unrefactored source.

