#!/usr/bin/env Rscript
# 01_pca_expansion.R ----
# PCA expansion for GSE118553 feature space.
#
# Inputs (objects expected to exist in the environment BEFORE sourcing/running):
#   - expr_mapped : probe-level expression with rownames = probe IDs
#   - map_df      : data.frame with at least columns: probe_id, symbol (gene symbol)
#   - expr_gene   : *core* gene expression matrix (genes x samples) used to define PCA sample order
#   - core_features : character vector of core gene symbols
#
# Outputs (written to outputs/ by default):
#   - expanded_genes.txt
#   - expanded_expression_gene_x_sample.csv
#   - pca_expansion_diagnostics.rds

suppressPackageStartupMessages({
  library(limma)
})

source(file.path("R", "utils","PCA_WGCNA.R"))

# Required packages
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

# ------------------------------------------------------------------
# Assumptions:
# The following objects must already exist in the environment:
#
# expr_mapped : probe x sample expression matrix (rownames = probe IDs)
# map_df      : data.frame with columns: probe, symbol
#               (symbol should be cleaned gene symbols)
# deg1, deg2, deg3 : limma topTable outputs
#               (deg*$feature = probe ID)
# ------------------------------------------------------------------

# 1) Define core probes based on overlap across contrasts
core_probes <- intersect(deg1$feature, union(deg2$feature, deg3$feature))

# 2) Convert probe IDs to gene symbols (core_features must be gene-level)
core_features <- map_df$symbol[match(core_probes, map_df$probe)]
core_features <- unique(core_features[!is.na(core_features) & core_features != ""])

cat("[core] probes:", length(core_probes), "\n")
cat("[core] gene symbols:", length(core_features), "\n")

# 3) Collapse probe-level expression to gene-level expression
#    If multiple probes map to the same gene, average expression values
expr_gene_all <- limma::avereps(expr_mapped, ID = map_df$symbol)

# 4) Subset gene-level expression matrix to core genes only
#    Final format: genes x samples matrix
core_features_in <- intersect(core_features, rownames(expr_gene_all))
expr_gene <- expr_gene_all[core_features_in, , drop = FALSE]

cat("[expr_gene] dim (genes x samples):", 
    paste(dim(expr_gene), collapse = " x "), "\n")

# (Optional) Add compatibility column if downstream scripts expect `probe_id`
if (!("probe_id" %in% colnames(map_df))) {
  map_df$probe_id <- map_df$probe
}


log_step("PCA expansion | input check")
require_objects(c("expr_mapped", "map_df", "expr_gene", "core_features"))

if (!("symbol" %in% colnames(map_df))) stop("map_df must contain a column named 'symbol' (gene symbol).")
if (is.null(rownames(expr_mapped))) stop("expr_mapped must have rownames = probe IDs.")

# --- 0) Make gene-level matrix from probe-level
log_step("PCA expansion | probe -> gene collapse (avereps)")

if (!is.matrix(expr_gene)) {
  # If expr_gene is a data.frame/tibble with a gene column, fix it
  if ("GeneSymbol" %in% colnames(expr_gene)) {
    rownames(expr_gene) <- expr_gene$GeneSymbol
    expr_gene$GeneSymbol <- NULL
  }
  expr_gene <- as.matrix(expr_gene)
}

# --- Align samples safely using intersection
common_samples <- intersect(colnames(expr_gene_all), colnames(expr_gene))

if (length(common_samples) < 5) {
  stop("Too few overlapping samples between expr_gene_all and expr_gene. Check sample IDs/colnames.")
}

# Reorder both matrices to the same sample order
expr_gene_all <- expr_gene_all[, common_samples, drop = FALSE]
expr_gene     <- expr_gene[,     common_samples, drop = FALSE]

log_step(sprintf("Aligned samples: %d common samples", length(common_samples)))

# Align samples to the core matrix sample order (very important)
expr_gene_all <- expr_gene_all[, colnames(expr_gene), drop = FALSE]

# --- 1) PCA on core genes
log_step("PCA expansion | PCA on core genes")

core_features <- unique(core_features)
core_features <- intersect(core_features, rownames(expr_gene_all))
if (length(core_features) < 10) stop("Too few core_features found in expr_gene_all after mapping (<10).")

X_core <- t(expr_gene_all[core_features, , drop = FALSE])  # samples x core genes
X_core <- safe_scale(X_core)

pca <- prcomp(X_core, center = TRUE, scale. = TRUE)

# --- Choose number of PCs by cumulative variance (statistically grounded)
K_var <- summary(pca)$importance[2, ]  # proportion variance
K_fixed <- which(cumsum(K_var) >= 0.80)[1]
K_fixed <- max(6, min(K_fixed, 25))

drop_PC1 <- FALSE 

use_idx <- seq_len(min(K_fixed, ncol(pca$x)))
if (drop_PC1) use_idx <- setdiff(use_idx, 1)

scores <- pca$x[, use_idx, drop = FALSE]  # samples x K_used
cat("Using PCs:", paste(colnames(scores), collapse = ", "), "\n")

# --- define m_and AFTER scores exists
m_and <- 2
alpha <- 0.05
alpha_soft <- 0.10
target_max <- 3000 

# --- 2) Correlations
log_step("PCA expansion | correlation with PCs (simple + partial)")

expr_gene_all <- expr_gene_all[, rownames(scores), drop = FALSE]  # align samples
X_all <- safe_scale(t(expr_gene_all))  # samples x genes

cor_with_pc <- function(x, y) suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))

# Simple correlations: genes x K
Cmat <- sapply(seq_len(ncol(scores)), function(k) {
  apply(X_all, 2, function(g) cor_with_pc(g, scores[, k]))
})
Cmat <- as.matrix(Cmat)
colnames(Cmat) <- colnames(scores)
rownames(Cmat) <- colnames(X_all)

n <- nrow(scores)

cor_to_p <- function(r, n) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  t <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t), df = n - 2)
}

P_simple <- apply(Cmat, 2, function(r) cor_to_p(r, n))
P_simple <- as.matrix(P_simple)
rownames(P_simple) <- rownames(Cmat)
colnames(P_simple) <- colnames(Cmat)

FDR_simple <- apply(P_simple, 2, p.adjust, method = "BH")
FDR_simple <- as.matrix(FDR_simple)
rownames(FDR_simple) <- rownames(Cmat)
colnames(FDR_simple) <- colnames(Cmat)

# Partial correlations (residual approach)
partial_cor_onegene <- function(g, S) {
  K <- ncol(S)
  out <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    if (K == 1) {
      out[k] <- cor_with_pc(g, S[, k])
    } else {
      others <- S[, -k, drop = FALSE]
      rg <- residuals(lm(g ~ others))
      rk <- residuals(lm(S[, k] ~ others))
      out[k] <- cor_with_pc(rg, rk)
    }
  }
  out
}

PCmat <- t(apply(X_all, 2, partial_cor_onegene, S = scores))
colnames(PCmat) <- colnames(scores)
rownames(PCmat) <- colnames(X_all)

df_partial <- n - (ncol(scores) - 1) - 2
df_partial <- max(df_partial, 5)

P_partial <- apply(PCmat, 2, function(r) cor_to_p(r, df_partial + 2))
P_partial <- as.matrix(P_partial)
rownames(P_partial) <- rownames(PCmat)
colnames(P_partial) <- colnames(PCmat)

FDR_partial <- apply(P_partial, 2, p.adjust, method = "BH")
FDR_partial <- as.matrix(FDR_partial)
rownames(FDR_partial) <- rownames(PCmat)
colnames(FDR_partial) <- colnames(PCmat)

# --- 3) FDR-based expansion (statistically significant genes)
log_step("PCA expansion | FDR-based expansion")

pass_p <- rowSums(FDR_partial < alpha, na.rm = TRUE) >= m_and
pass_s <- rowSums(FDR_simple  < alpha_soft, na.rm = TRUE) >= 1  # simple은 보조 확인

keep <- pass_p & pass_s

cand_genes <- rownames(Cmat)[keep]
cand_genes <- setdiff(cand_genes, core_features)

expanded_genes <- unique(c(core_features, cand_genes))

# Optional: only CAP (do not fill). If expanded is too large, cut by combined evidence score.
if (length(expanded_genes) > target_max) {
  pool <- setdiff(expanded_genes, core_features)
  
  score_partial <- apply(abs(PCmat[pool, , drop = FALSE]), 1, function(v)
    mean(sort(v, decreasing = TRUE)[1:m_and])
  )
  score_simple <- apply(abs(Cmat[pool, , drop = FALSE]), 1, function(v)
    mean(sort(v, decreasing = TRUE)[1:m_and])
  )
  score <- score_partial * score_simple
  
  ranked <- names(sort(score, decreasing = TRUE))
  keep_add <- head(ranked, target_max - length(core_features))
  expanded_genes <- unique(c(core_features, keep_add))
}

expanded_expr  <- expr_gene_all[expanded_genes, , drop = FALSE]

cat("alpha =", alpha, "alpha_soft =", alpha_soft, "m_and =", m_and, "\n")
cat("Cand genes (non-core):", length(cand_genes), "\n")
cat("FINAL expanded genes:", length(expanded_genes), "\n")
cat("expanded_expr dim:", nrow(expanded_expr), "x", ncol(expanded_expr), "\n")

best <- list(
  method = "FDR_expansion",
  K_fixed = K_fixed,
  drop_PC1 = drop_PC1,
  m_and = m_and,
  alpha = alpha,
  alpha_soft = alpha_soft,
  expanded_n = length(expanded_genes)
)

cat("Chosen thresholds: thr_p =", best$thr_p, "thr_s =", best$thr_s, "\n")
cat("Genes after AND filter (including core union):", best$nn, "\n")

cat("FINAL expanded genes:", length(expanded_genes), "\n")
cat("expanded_expr dim:", nrow(expanded_expr), "x", ncol(expanded_expr), "\n")

# --- 5) Save outputs
log_step("PCA expansion | save outputs")
out_dir <- "outputs"
write_lines(expanded_genes, file.path(out_dir, "expanded_genes.txt"))
write_csv(
  data.frame(GeneSymbol = rownames(expanded_expr), expanded_expr, check.names = FALSE),
  file.path(out_dir, "expanded_expression_gene_x_sample.csv")
)

saveRDS(
  list(
    params = list(K_fixed = K_fixed, drop_PC1 = drop_PC1, m_and = m_and, target_max = target_max,
                  thr_p = best$thr_p, thr_s = best$thr_s),
    pca = pca,
    scores = scores,
    Cmat = Cmat,
    PCmat = PCmat
  ),
  file = file.path(out_dir, "pca_expansion_diagnostics.rds")
)

log_step("PCA expansion | done")
