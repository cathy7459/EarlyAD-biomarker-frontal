log_msg("RUNNING FILE: scripts/03b_enet_train_validate.R")

# ============================================================
# Paths
# ============================================================
deg_path <- file.path("outputs", "tables", "deg_gse33000_postproj.csv")
deg_f_path <- file.path("outputs", "tables", "deg_filtered_genes_postproj.csv")
expr_path <- file.path("data", "processed", "gse33000_expr_filtered_gene_x_sample.csv")

if (!file.exists(deg_path)) {
  stop("Missing DEG table: ", deg_path,
       "\nRun scripts/02_deg_effect_filter.R first.")
}
if (!file.exists(deg_f_path)) {
  stop("Missing filtered DEG table: ", deg_f_path,
       "\nRun scripts/02_deg_effect_filter.R first.")
}
if (!file.exists(expr_path)) {
  stop("Missing filtered expression matrix: ", expr_path,
       "\nRun scripts/02_deg_effect_filter.R first.")
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
# 1) Load DEG tables + filtered expression matrix
# ============================================================
deg_dt <- data.table::fread(deg_path)
deg_f_dt <- data.table::fread(deg_f_path)
expr_dt <- data.table::fread(expr_path)

stopifnot("GeneSymbol" %in% names(expr_dt))
stopifnot("GeneSymbol" %in% names(deg_dt))
stopifnot("GeneSymbol" %in% names(deg_f_dt))

expr <- as.data.frame(expr_dt[, -1, with = FALSE])
rownames(expr) <- trimws(as.character(expr_dt$GeneSymbol))
expr[] <- lapply(expr, as.numeric)

deg_keep <- unique(trimws(as.character(deg_f_dt$GeneSymbol)))
deg_keep <- deg_keep[!is.na(deg_keep) & deg_keep != ""]

common_deg_genes <- intersect(rownames(expr), deg_keep)
if (length(common_deg_genes) < 5) {
  stop("Too few overlap genes between expr_path and deg_f_path: ", length(common_deg_genes))
}
expr <- expr[common_deg_genes, , drop = FALSE]

log_msg("[1/7] Loaded DEG table: ", nrow(deg_dt))
log_msg("[1/7] Loaded filtered DEG table: ", nrow(deg_f_dt))
log_msg("[1/7] Training expression genes after DEG filter: ", nrow(expr))

expr_sample_ids_raw <- colnames(expr)
expr_sample_ids_norm <- normalize_sample_ids(expr_sample_ids_raw)
keep_expr_cols <- !duplicated(expr_sample_ids_norm)
expr <- expr[, keep_expr_cols, drop = FALSE]
expr_sample_ids_raw <- expr_sample_ids_raw[keep_expr_cols]
expr_sample_ids_norm <- expr_sample_ids_norm[keep_expr_cols]
colnames(expr) <- expr_sample_ids_norm

core_genes <- NULL
if (!is.null(state$core_genes) && length(state$core_genes) > 0) {
  core_genes <- unique(trimws(as.character(state$core_genes)))
  core_genes <- core_genes[!is.na(core_genes) & core_genes != ""]
}
if (!is.null(core_genes)) {
  common_core <- intersect(rownames(expr), core_genes)
  if (length(common_core) >= 20) {
    expr <- expr[common_core, , drop = FALSE]
    log_msg("[1/7] Restricted to core genes: ", nrow(expr))
  }
}

# ============================================================
# 2) Load metadata and align
# ============================================================
meta <- as.data.frame(data.table::fread(meta_path))
sample_id_candidates <- c("sample_id", "geo_accession", "geo accession", "SampleID", "sample", "gsm", "GSM", "Sample")
sample_id_col <- sample_id_candidates[sample_id_candidates %in% names(meta)][1]
if (is.na(sample_id_col) || !nzchar(sample_id_col)) {
  stop("Could not find sample ID column in metadata.\nAvailable columns:\n", paste(names(meta), collapse = ", "))
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
x_raw <- t(as.matrix(expr))
storage.mode(x_raw) <- "double"
rownames(x_raw) <- colnames(expr)

if (anyNA(x_raw)) {
  log_msg("[3/7] NA detected -> median imputation")
  for (j in seq_len(ncol(x_raw))) {
    idx <- is.na(x_raw[, j])
    if (any(idx)) x_raw[idx, j] <- median(x_raw[, j], na.rm = TRUE)
  }
}

y <- factor(meta$condition, levels = c("Control", "AD"))
if (anyNA(y)) stop("Outcome contains NA values.")
y_bin <- ifelse(y == "AD", 1, 0)

# Remove non-constant features using raw training data
feature_sd <- apply(x_raw, 2, sd, na.rm = TRUE)
keep_feat <- is.finite(feature_sd) & !is.na(feature_sd) & feature_sd > 0
x_raw <- x_raw[, keep_feat, drop = FALSE]
if (ncol(x_raw) < 2) stop("Too few non-constant features after filtering.")

# Save training-side statistics on the same feature space used for modeling
feature_medians_full <- apply(x_raw, 2, median, na.rm = TRUE)
feature_means_full   <- colMeans(x_raw, na.rm = TRUE)
feature_sds_full     <- apply(x_raw, 2, sd, na.rm = TRUE)

feature_medians_full[!is.finite(feature_medians_full)] <- 0
feature_means_full[!is.finite(feature_means_full)] <- 0
feature_sds_full[!is.finite(feature_sds_full) | feature_sds_full == 0] <- 1

# IMPORTANT:
# Standardize training data using training means and SDs,
# so that training and external test preprocessing are identical.
x <- sweep(x_raw, 2, feature_means_full, FUN = "-")
x <- sweep(x, 2, feature_sds_full, FUN = "/")
storage.mode(x) <- "double"

if (any(!is.finite(x))) {
  bad_n <- sum(!is.finite(x))
  log_msg("[3/7] Non-finite values detected after training scaling: ", bad_n, " -> replacing with 0")
  x[!is.finite(x)] <- 0
}

log_msg("[3/7] Input ready: samples=", nrow(x), ", genes=", ncol(x), " (training data standardized)")

# ------------------------------------------------------------
# Optional: prune highly correlated hub/core genes
# ------------------------------------------------------------
prune_correlated_features <- function(x, cor_cutoff = 0.90) {
  x <- as.matrix(x)
  cmat <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
  cmat[is.na(cmat)] <- 0
  
  keep <- rep(TRUE, ncol(cmat))
  names(keep) <- colnames(cmat)
  
  for (i in seq_len(ncol(cmat) - 1)) {
    if (!keep[i]) next
    high_cor_idx <- which(abs(cmat[i, ]) > cor_cutoff)
    high_cor_idx <- high_cor_idx[high_cor_idx > i]
    if (length(high_cor_idx) > 0) {
      keep[high_cor_idx] <- FALSE
    }
  }
  
  colnames(x)[keep]
}

# Apply only when a core/hub gene set is being used
if (!is.null(core_genes)) {
  genes_before_prune <- colnames(x)
  pruned_genes <- prune_correlated_features(x[, intersect(colnames(x), core_genes), drop = FALSE], cor_cutoff = 0.90)
  
  if (length(pruned_genes) >= 10) {
    x <- x[, pruned_genes, drop = FALSE]
    log_msg("[3/7] Correlation-pruned hub/core genes: ", ncol(x))
  }
}

# ============================================================
# 4) Univariate prefilter
# ============================================================
univ_res <- run_univariate_filter(
  x = x,
  y_bin = y_bin,
  method = UNIV_METHOD,
  use_fdr_gate = USE_FDR_GATE,
  fdr_cutoff = FDR_CUTOFF,
  min_fdr_pass = MIN_FDR_PASS,
  top_n = UNIV_TOP_N
)
univ_dt <- univ_res$ranking
prefilter_genes <- univ_res$selected_features
if (length(prefilter_genes) < 5) stop("Too few genes selected by univariate prefilter.")
x_filt <- x[, prefilter_genes, drop = FALSE]
data.table::fwrite(univ_dt, out_univ)
log_msg("[4/7] Prefiltered genes: ", ncol(x_filt))

# ============================================================
# 5) Repeated cross-validation for Elastic Net alpha selection
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
# 6) Final ENet fit / selected genes / predictions
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

# ============================================================
# 7) Save model bundle
# ============================================================
bundle <- list(
  created_at = Sys.time(),
  script = "03b_enet_train_validate.R",
  cvfit = cvfit,
  deg_table_path = deg_path,
  deg_filtered_path = deg_f_path,
  expr_filtered_path = expr_path,
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
