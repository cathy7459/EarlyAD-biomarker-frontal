#!/usr/bin/env Rscript
# run_all.R ----
# Minimal runner. This script assumes you already ran your upstream preprocessing
# to create objects: expr_mapped, map_df, expr_gene (core), core_features, (optional) group_f.
#
# Typical usage:
#   1) Source your preprocessing code (that creates the required objects)
#   2) source("scripts/run_all.R")  OR  Rscript scripts/run_all.R after saving .RData

source(file.path("R/utils", PCA_WGCNA.R"))

log_step("Run all | note")
cat("This repo starts from objects produced upstream (preprocessing / GEO import).\n",
    "Required: expr_mapped, map_df, expr_gene (core), core_features.\n",
    "Optional: group_f (named by sample IDs).\n", sep = "")

# PCA expansion
source(file.path("scripts", "01_pca_expansion.R"))

# After PCA expansion, load expanded matrix back as expr_gene for WGCNA
expanded_path <- file.path("outputs", "expanded_expression_gene_x_sample.csv")
if (!file.exists(expanded_path)) stop("Expected output not found: ", expanded_path)

log_step("Run all | load expanded matrix for WGCNA")
expanded_df <- read.csv(expanded_path, check.names = FALSE)
rownames(expanded_df) <- expanded_df$GeneSymbol
expanded_df$GeneSymbol <- NULL
expr_gene <- as.matrix(expanded_df)

# WGCNA
source(file.path("scripts", "02_wgcna.R"))

log_step("Run all | finished")
