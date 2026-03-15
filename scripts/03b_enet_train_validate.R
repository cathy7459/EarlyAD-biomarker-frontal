log_msg("RUNNING FILE: scripts/03b_enet_train_validate.R")

# ============================================================
# 0) Paths
# ============================================================
expr_path <- file.path("data", "processed", "gse33000_expr_projected_gene_x_sample.csv")
if (!file.exists(expr_path)) {
  stop(
    "Missing projected training expression matrix: ", expr_path,
    "\nRun scripts/01_build_expr_gse33000_with_gpl.R first."
  )
}

meta_path <- file.path("data", "processed", "GSE33000_pheno.csv")
if (!file.exists(meta_path)) stop("Missing metadata: ", meta_path)

out_tables <- file.path(CFG$out_dir, "tables")
out_figures <- file.path(CFG$out_dir, "figures")
.dir_create(out_tables); .dir_create(out_figures)

bundle_prefix <- if (!is.null(CFG$prefix) && nzchar(CFG$prefix)) CFG$prefix else "ad_enet"
out_sel        <- file.path(out_tables, paste0(bundle_prefix, "_enet_selected_genes.csv"))
out_coef       <- file.path(out_tables, paste0(bundle_prefix, "_enet_coefficients.csv"))
out_pred_oof   <- file.path(out_tables, paste0(bundle_prefix, "_enet_oof_predictions.csv"))
out_pred_app   <- file.path(out_tables, paste0(bundle_prefix, "_enet_apparent_predictions.csv"))
out_alpha      <- file.path(out_tables, paste0(bundle_prefix, "_enet_alpha_repeatcv_summary.csv"))
out_alpha_raw  <- file.path(out_tables, paste0(bundle_prefix, "_enet_alpha_repeatcv_raw.csv"))
out_univ       <- file.path(out_tables, paste0(bundle_prefix, "_univariate_filter_ranking.csv"))
out_roc_oof    <- file.path(out_figures, paste0(bundle_prefix, "_roc_train_oof_enet.png"))
out_roc_app    <- file.path(out_figures, paste0(bundle_prefix, "_roc_train_apparent_enet.png"))
out_bundle     <- CFG$model_bundle_rds
.dir_create(dirname(out_bundle))

# ============================================================
# 1) Load projected expression matrix
# ============================================================
expr_dt <- data.table::fread(expr_path)
stopifnot("GeneSymbol" %in% names(expr_dt))

expr <- as.data.frame(expr_dt[, -1, with = FALSE])
rownames(expr) <- trimws(as.character(expr_dt$GeneSymbol))
expr[] <- lapply(expr, as.numeric)

log_msg("[1/7] Loaded projected training matrix: ", nrow(expr), " genes")

# ============================================================
# 2) Load metadata and align samples
# ============================================================
meta <- as.data.frame(data.table::fread(meta_path))

sample_id_candidates <- c(
  "sample_id", "geo_accession", "geo accession",
  "SampleID", "sample", "gsm", "GSM", "Sample"
)
sample_id_col <- sample_id_candidates[sample_id_candidates %in% names(meta)][1]

if (is.na(sample_id_col) || !nzchar(sample_id_col)) {
  stop(
    "Could not find sample ID column in metadata.\nAvailable columns:\n",
    paste(names(meta), collapse = ", ")
  )
}

if (!"condition" %in% names(meta)) {
  cond_raw <- apply(meta, 1, function(x) paste(as.character(x), collapse = " | "))
  meta$condition <- ifelse(
    grepl("alzheimer's disease|alzheimer", cond_raw, ignore.case = TRUE), "AD",
    ifelse(grepl("non-demented|control", cond_raw, ignore.case = TRUE), "Control", NA_character_)
  )
}

meta$sample_id_raw  <- trimws(as.character(meta[[sample_id_col]]))
meta$sample_id_norm <- normalize_sample_ids(meta$sample_id_raw)
meta$condition      <- trimws(as.character(meta$condition))

meta <- meta[meta$condition %in% c("AD", "Control"), , drop = FALSE]
meta <- meta[!duplicated(meta$sample_id_norm), , drop = FALSE]

expr_sample_ids_raw  <- colnames(expr)
expr_sample_ids_norm <- normalize_sample_ids(expr_sample_ids_raw)

keep_expr_cols <- !duplicated(expr_sample_ids_norm)
expr <- expr[, keep_expr_cols, drop = FALSE]
expr_sample_ids_raw  <- expr_sample_ids_raw[keep_expr_cols]
expr_sample_ids_norm <- expr_sample_ids_norm[keep_expr_cols]
colnames(expr) <- expr_sample_ids_norm

common_samples <- intersect(colnames(expr), meta$sample_id_norm)
if (length(common_samples) < 20) {
  stop("Too few matched samples between expression and metadata: ", length(common_samples))
}

expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id_norm), , drop = FALSE]
stopifnot(identical(colnames(expr), meta$sample_id_norm))

colnames(expr) <- meta$sample_id_raw

log_msg("[2/7] Matched samples: ", ncol(expr))
log_msg("[2/7] AD samples: ", sum(meta$condition == "AD"))
log_msg("[2/7] Control samples: ", sum(meta$condition == "Control"))

# ============================================================
# 3) Build model input
# ============================================================
x <- t(as.matrix(expr))
storage.mode(x) <- "double"
rownames(x) <- colnames(expr)

# Remove unusable features before imputation
non_all_na <- colSums(!is.na(x)) > 0
x <- x[, non_all_na, drop = FALSE]

gene_sds0 <- apply(x, 2, sd, na.rm = TRUE)
keep_var0 <- !is.na(gene_sds0) & gene_sds0 > 0
x <- x[, keep_var0, drop = FALSE]

log_msg("[3/7] Features after NA/variance filtering: ", ncol(x))

if (anyNA(x)) {
  log_msg("[3/7] Remaining NA detected -> column-median imputation")
  na_cols <- which(colSums(is.na(x)) > 0)
  col_medians <- apply(x[, na_cols, drop = FALSE], 2, median, na.rm = TRUE)
  for (jj in seq_along(na_cols)) {
    j <- na_cols[jj]
    idx <- is.na(x[, j])
    if (any(idx)) {
      x[idx, j] <- col_medians[jj]
    }
  }
} else {
  log_msg("[3/7] No missing values remained after filtering")
}

y <- factor(meta$condition, levels = c("Control", "AD"))
y_bin <- ifelse(y == "AD", 1, 0)

# Standardize features before correlation pruning and Elastic Net
feature_means <- colMeans(x, na.rm = TRUE)
feature_sds   <- apply(x, 2, sd, na.rm = TRUE)
feature_sds[is.na(feature_sds) | feature_sds == 0] <- 1

x_scaled <- sweep(x, 2, feature_means, FUN = "-")
x_scaled <- sweep(x_scaled, 2, feature_sds, FUN = "/")

# ============================================================
# 4) Correlation pruning
# ============================================================
pruned_genes <- prune_correlated_features(x_scaled, cor_cutoff = 0.90)

if (length(pruned_genes) < 10) {
  log_msg("[4/7] Correlation pruning retained too few genes -> using unpruned projected features")
  x_pruned <- x_scaled
} else {
  x_pruned <- x_scaled[, pruned_genes, drop = FALSE]
  log_msg("[4/7] Correlation-pruned projected features: ", ncol(x_pruned))
}

# ============================================================
# 5) Univariate prefilter
# ============================================================
univ_res <- run_univariate_filter(
  x = x_pruned,
  y_bin = y_bin,
  method = UNIV_METHOD,
  use_fdr_gate = USE_FDR_GATE,
  fdr_cutoff = FDR_CUTOFF,
  min_fdr_pass = MIN_FDR_PASS,
  top_n = UNIV_TOP_N
)

univ_dt <- univ_res$ranking
prefilter_genes <- univ_res$selected_features

data.table::fwrite(univ_dt, out_univ)

if (length(prefilter_genes) < 5) {
  stop("Too few genes remained after univariate prefilter: ", length(prefilter_genes))
}

x_filt <- x_pruned[, prefilter_genes, drop = FALSE]

log_msg("[5/7] Univariate prefilter retained: ", ncol(x_filt), " genes")


# ============================================================
# 6) Repeated cross-validation for Elastic Net alpha selection
# ============================================================

# Evaluate model performance across multiple CV repeats
# to obtain a stable estimate of out-of-fold AUC.

perf_list <- vector("list", length(ALPHA_GRID) * CV_REPEATS)
kk <- 1

for (rr in seq_len(CV_REPEATS)) {
  
  # Generate stratified fold assignments
  foldid <- make_stratified_foldid(y_bin, k = CV_FOLDS, seed = SEED + rr)
  
  for (aa in ALPHA_GRID) {
    
    # Fit Elastic Net model with the given alpha
    cvfit_rr <- fit_cv_glmnet(x_filt, y_bin, alpha = aa, foldid = foldid)
    
    # Select regularization parameter using lambda.1se
    lam_rr <- get_lambda_1se(cvfit_rr)
    
    # Obtain out-of-fold predicted probabilities
    prob_oof_rr <- get_oof_predictions(
      x_filt,
      y_bin,
      alpha = aa,
      foldid = foldid,
      lambda_value = lam_rr
    )
    
    # Compute out-of-fold AUC
    auc_oof_rr <- as.numeric(
      pROC::auc(pROC::roc(y_bin, prob_oof_rr, quiet = TRUE))
    )
    
    # Extract selected genes for this model
    coef_rr <- extract_coef_table(cvfit_rr, lam_rr)
    sel_rr  <- extract_selected_genes(coef_rr)
    
    perf_list[[kk]] <- data.table::data.table(
      repeat_id   = rr,
      alpha       = aa,
      lambda_1se  = lam_rr,
      oof_auc     = auc_oof_rr,
      n_genes     = nrow(sel_rr)
    )
    
    kk <- kk + 1
  }
}

perf_dt <- data.table::rbindlist(perf_list)

# Save raw cross-validation results
data.table::fwrite(perf_dt, out_alpha_raw)


# ============================================================
# Summarize performance across repeats
# ============================================================

alpha_summary <- perf_dt[, .(
  mean_oof_auc   = mean(oof_auc, na.rm = TRUE),
  sd_oof_auc     = sd(oof_auc, na.rm = TRUE),
  mean_n_genes   = mean(n_genes, na.rm = TRUE),
  median_n_genes = median(n_genes, na.rm = TRUE)
), by = alpha]


# ============================================================
# Alpha selection strategy
# ============================================================

# Step 1
# Identify models with performance close to the best AUC.
# A small tolerance is allowed to avoid over-selecting LASSO (alpha = 1).

best_auc <- max(alpha_summary$mean_oof_auc, na.rm = TRUE)

ALPHA_TOL <- 0.01

alpha_keep <- alpha_summary[
  mean_oof_auc >= (best_auc - ALPHA_TOL)
]

if (nrow(alpha_keep) == 0) {
  alpha_keep <- alpha_summary[which.max(mean_oof_auc)]
}


# Step 2
# Apply a mild preference toward intermediate alpha values
# to encourage true Elastic Net behavior rather than pure LASSO.

TARGET_ALPHA <- 0.70
ALPHA_CENTER_PENALTY <- 0.002


# Step 3
# Apply a weak complexity penalty based on the number of selected genes.
# This prevents extremely large models while avoiding strong bias toward
# overly sparse models.

ALPHA_COMPLEXITY_PENALTY <- 0.0005


alpha_keep[, score :=
             mean_oof_auc
           - ALPHA_COMPLEXITY_PENALTY * mean_n_genes
           - ALPHA_CENTER_PENALTY * abs(alpha - TARGET_ALPHA)
]


# Rank candidate models using the combined score
data.table::setorder(alpha_keep, -score, -mean_oof_auc, mean_n_genes)


# Final alpha selection
best_alpha <- alpha_keep$alpha[1]


# Record selected alpha
alpha_summary[, selected := alpha == best_alpha]

alpha_summary[, selection_score := NA_real_]

alpha_summary[
  alpha_keep,
  selection_score := i.score,
  on = "alpha"
]

data.table::setorder(alpha_summary, -mean_oof_auc, mean_n_genes)

# Save alpha selection summary
data.table::fwrite(alpha_summary, out_alpha)

log_msg("[5/7] Selected alpha: ", best_alpha)

# ============================================================
# 7) Final ENet fit / selected genes / predictions
# ============================================================
foldid_final <- make_stratified_foldid(y_bin, k = CV_FOLDS, seed = SEED)
cvfit <- fit_cv_glmnet(x_filt, y_bin, alpha = best_alpha, foldid = foldid_final)
lam <- get_lambda_1se(cvfit)
coef_dt <- extract_coef_table(cvfit, lam)
sel_dt <- extract_selected_genes(coef_dt)

if (nrow(sel_dt) > MAX_FINAL_GENES) {
  sel_keep <- sel_dt$Feature[seq_len(MAX_FINAL_GENES)]
  x_filt2 <- x_filt[, sel_keep, drop = FALSE]
  cvfit <- fit_cv_glmnet(x_filt2, y_bin, alpha = best_alpha, foldid = foldid_final)
  lam <- get_lambda_1se(cvfit)
  coef_dt <- extract_coef_table(cvfit, lam)
  sel_dt <- extract_selected_genes(coef_dt)
  x_filt <- x_filt2
}

if (nrow(sel_dt) < MIN_FINAL_GENES) {
  warning("Final ENet selected fewer genes than MIN_FINAL_GENES.")
}

data.table::fwrite(sel_dt, out_sel)
data.table::fwrite(coef_dt, out_coef)

prob_app <- as.numeric(predict(cvfit, newx = x_filt, s = lam, type = "response"))
prob_oof <- get_oof_predictions(x_filt, y_bin, alpha = best_alpha, foldid = foldid_final, lambda_value = lam)
auc_app <- as.numeric(pROC::auc(pROC::roc(y_bin, prob_app, quiet = TRUE)))
auc_oof <- as.numeric(pROC::auc(pROC::roc(y_bin, prob_oof, quiet = TRUE)))

pred_app_dt <- data.table::data.table(sample_id = rownames(x_filt), condition = as.character(y), y_bin = y_bin, prob_AD = prob_app)
pred_oof_dt <- data.table::data.table(sample_id = rownames(x_filt), condition = as.character(y), y_bin = y_bin, prob_AD = prob_oof)
data.table::fwrite(pred_app_dt, out_pred_app)
data.table::fwrite(pred_oof_dt, out_pred_oof)

roc_app <- pROC::roc(y_bin, prob_app, quiet = TRUE)
roc_oof <- pROC::roc(y_bin, prob_oof, quiet = TRUE)
png(out_roc_app, width = 1800, height = 1600, res = 220)
plot(roc_app, main = sprintf("Training Apparent ROC (AUC = %.3f)", auc_app), legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 2)
dev.off()
png(out_roc_oof, width = 1800, height = 1600, res = 220)
plot(roc_oof, main = sprintf("Training ROC (AUC = %.3f)", auc_oof), legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 2)
dev.off()

feature_medians_full <- apply(x_filt, 2, median, na.rm = TRUE)
feature_means_full   <- colMeans(x_filt, na.rm = TRUE)
feature_sds_full     <- apply(x_filt, 2, sd, na.rm = TRUE)
feature_sds_full[is.na(feature_sds_full) | feature_sds_full == 0] <- 1
# ============================================================
# 8) Save model bundle
# ============================================================
bundle <- list(
  created_at = Sys.time(),
  script = "03b_enet_train_validate.R",
  cvfit = cvfit,
  expr_train_path = expr_path,
  lambda_use = lam,
  lambda_mode = LAMBDA_MODE,
  alpha_use = best_alpha,
  genes_used = colnames(x_filt),
  prefilter_genes = prefilter_genes,
  selected = sel_dt,
  coef_table = coef_dt,
  alpha_summary = alpha_summary,
  alpha_repeat_perf = perf_dt,
  univariate_ranking = univ_dt,
  feature_medians = feature_medians_full,
  feature_means = feature_means_full,
  feature_sds = feature_sds_full,
  n_samples = nrow(x_filt),
  n_genes_model = ncol(x_filt),
  sample_ids = rownames(x_filt),
  y_train = y_bin,
  y_train_chr = as.character(y),
  apparent_pred = prob_app,
  oof_pred = prob_oof,
  apparent_auc = auc_app,
  oof_auc = auc_oof,
  seed = SEED,
  cv_folds = CV_FOLDS,
  cv_repeats = CV_REPEATS,
  univ_top_n = UNIV_TOP_N,
  univ_method = UNIV_METHOD,
  use_fdr_gate = USE_FDR_GATE,
  fdr_cutoff = FDR_CUTOFF,
  min_fdr_pass = MIN_FDR_PASS
)

saveRDS(bundle, out_bundle)
state$train_x_filt <- x_filt
state$train_y_bin <- y_bin
state$train_meta <- meta
state$enet_bundle <- bundle

log_msg("DONE. samples=", nrow(x_filt),
        " genes=", ncol(x_filt),
        " selected=", nrow(sel_dt),
        " apparent_AUC=", round(auc_app, 4),
        " OOF_AUC=", round(auc_oof, 4),
        " bundle=", out_bundle)
