# EarlyAD-biomarker-PCA-WGCNA
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
