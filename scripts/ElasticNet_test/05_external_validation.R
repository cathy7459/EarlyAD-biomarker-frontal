# scripts/05_external_validate_gse122063_frontal.R

log_msg("RUNNING FILE: scripts/05_external_validate_gse122063_frontal.R")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

if (!exists("state")) state <- new.env(parent = emptyenv())
if (!exists("CFG")) stop("CFG not found. source('config/config.R') first.")

# ============================================================
# 0) Settings
# ============================================================
TEST_REGION_USE <- "frontal"

out_tables  <- file.path(CFG$out_dir, "tables")
out_figures <- file.path(CFG$out_dir, "figures")
dir.create(out_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)

bundle_prefix <- if (!is.null(CFG$prefix) && nzchar(CFG$prefix)) CFG$prefix else "ad_enet"
region_tag <- paste0(TEST_REGION_USE, "_excluding_vascular")

out_pred_ext <- file.path(
  out_tables,
  paste0(bundle_prefix, "_external_", CFG$test_gse_id, "_", region_tag, "_predictions.csv")
)

out_perf_ext <- file.path(
  out_tables,
  paste0(bundle_prefix, "_external_", CFG$test_gse_id, "_", region_tag, "_performance.csv")
)

out_roc_ext <- file.path(
  out_figures,
  paste0(bundle_prefix, "_external_", CFG$test_gse_id, "_", region_tag, "_roc.png")
)

# ============================================================
# 1) Load model bundle
# ============================================================
if (!file.exists(CFG$model_bundle_rds)) {
  stop("Model bundle not found: ", CFG$model_bundle_rds)
}

bundle <- readRDS(CFG$model_bundle_rds)

if (is.null(bundle$cvfit)) {
  stop("cvfit missing in model bundle.")
}
if (is.null(bundle$lambda_use) || !is.finite(bundle$lambda_use)) {
  stop("lambda_use missing or invalid in model bundle.")
}
if (is.null(bundle$genes_used) || length(bundle$genes_used) < 2) {
  stop("genes_used missing or too short in model bundle.")
}

model_fit <- bundle$cvfit$glmnet.fit
lambda_use <- bundle$lambda_use
genes_used <- as.character(bundle$genes_used)

log_msg("[1/5] Loaded model bundle: ", CFG$model_bundle_rds)
log_msg("[1/5] Model genes: ", length(genes_used))
log_msg("[1/5] Lambda used: ", signif(lambda_use, 6))
log_msg("[1/5] Alpha used: ", bundle$alpha_use)

# ============================================================
# 2) Load prepared external test set from state
#    and force re-align to current model feature space
# ============================================================
if (is.null(state$X_test)) {
  stop(
    "state$X_test not found.\n",
    "Run scripts/04_prepare_external_test_set.R first."
  )
}
if (is.null(state$y_test) || is.null(state$y_test_chr)) {
  stop(
    "state$y_test / state$y_test_chr not found.\n",
    "Run scripts/04_prepare_external_test_set.R first."
  )
}
if (is.null(state$test_region_used)) {
  stop(
    "state$test_region_used not found.\n",
    "Run scripts/04_prepare_external_test_set.R first."
  )
}

if (tolower(state$test_region_used) != TEST_REGION_USE) {
  stop(
    "Prepared test region is '", state$test_region_used,
    "' but this script expects '", TEST_REGION_USE, "'."
  )
}

X_test_raw <- as.matrix(state$X_test)
storage.mode(X_test_raw) <- "double"

y_test <- as.integer(state$y_test)
y_test_chr <- as.character(state$y_test_chr)

if (nrow(X_test_raw) != length(y_test)) {
  stop("Mismatch between nrow(X_test) and length(y_test).")
}
if (!all(y_test %in% c(0L, 1L))) {
  stop("y_test must contain only 0/1.")
}

if (is.null(colnames(X_test_raw)) || any(colnames(X_test_raw) == "")) {
  stop("state$X_test must have valid gene column names.")
}

# re-align test matrix to current bundle genes_used
X_test <- matrix(
  NA_real_,
  nrow = nrow(X_test_raw),
  ncol = length(genes_used),
  dimnames = list(rownames(X_test_raw), genes_used)
)

matched_cols <- intersect(colnames(X_test_raw), genes_used)
missing_cols <- setdiff(genes_used, colnames(X_test_raw))

if (length(matched_cols) > 0) {
  X_test[, matched_cols] <- X_test_raw[, matched_cols, drop = FALSE]
}

# impute missing columns using train medians from current bundle
if (length(missing_cols) > 0) {
  if (is.null(bundle$feature_medians)) {
    stop("bundle$feature_medians is missing; cannot impute missing test genes.")
  }
  
  fill_vals <- bundle$feature_medians[missing_cols]
  fill_vals[is.na(fill_vals) | !is.finite(fill_vals)] <- 0
  
  for (g in missing_cols) {
    X_test[, g] <- fill_vals[g]
  }
}

storage.mode(X_test) <- "double"

if (anyNA(X_test) || any(!is.finite(X_test))) {
  stop("X_test contains NA or non-finite values after re-alignment.")
}

missing_fraction <- length(missing_cols) / length(genes_used)

log_msg("[2/5] External test matrix re-aligned to current bundle")
log_msg("[2/5] matched genes: ", length(matched_cols), "/", length(genes_used))
log_msg("[2/5] missing genes imputed from train medians: ", length(missing_cols))
log_msg("[2/5] missing fraction: ", round(missing_fraction, 4))
log_msg("[2/5] External test matrix ready: ", nrow(X_test), " samples x ", ncol(X_test), " genes")
log_msg("[2/5] External labels: AD=", sum(y_test == 1), " | Control=", sum(y_test == 0))

if (missing_fraction > 0.2) {
  warning("More than 20% of model genes were absent in state$X_test and were median-imputed.")
}
# ============================================================
# 3) Predict on external test set
# ============================================================
prob_ext <- as.numeric(
  predict(
    model_fit,
    newx = X_test,
    s = lambda_use,
    type = "response"
  )
)

if (length(prob_ext) != nrow(X_test)) {
  stop("Prediction length mismatch.")
}
if (anyNA(prob_ext) || any(!is.finite(prob_ext))) {
  stop("External predictions contain NA or non-finite values.")
}

pred_class <- ifelse(prob_ext >= 0.5, 1L, 0L)

# ============================================================
# 4) Evaluate external performance
# ============================================================
roc_ext <- pROC::roc(
  response = y_test,
  predictor = prob_ext,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- as.numeric(pROC::auc(roc_ext))
ci_auc <- as.numeric(pROC::ci.auc(roc_ext))

coords_05 <- pROC::coords(
  roc_ext,
  x = 0.5,
  input = "threshold",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

acc_ext <- mean(pred_class == y_test)
sens_ext <- as.numeric(coords_05["sensitivity"])
spec_ext <- as.numeric(coords_05["specificity"])

tp <- sum(pred_class == 1 & y_test == 1)
tn <- sum(pred_class == 0 & y_test == 0)
fp <- sum(pred_class == 1 & y_test == 0)
fn <- sum(pred_class == 0 & y_test == 1)

pred_ext_dt <- data.table(
  sample_id = rownames(X_test),
  condition = y_test_chr,
  y_bin = y_test,
  prob_AD = prob_ext,
  pred_class_0.5 = pred_class
)

perf_ext_dt <- data.table(
  dataset = CFG$test_gse_id,
  region = TEST_REGION_USE,
  n_samples = nrow(X_test),
  n_AD = sum(y_test == 1),
  n_Control = sum(y_test == 0),
  auc = auc_ext,
  auc_ci_lower = ci_auc[1],
  auc_ci_mid = ci_auc[2],
  auc_ci_upper = ci_auc[3],
  accuracy_threshold_0_5 = acc_ext,
  sensitivity_threshold_0_5 = sens_ext,
  specificity_threshold_0_5 = spec_ext,
  TP = tp,
  TN = tn,
  FP = fp,
  FN = fn,
  model_alpha = bundle$alpha_use,
  model_lambda = lambda_use,
  n_model_genes = length(genes_used),
  n_selected_nonzero = if (!is.null(bundle$selected)) nrow(bundle$selected) else NA_integer_,
  train_apparent_auc = if (!is.null(bundle$apparent_auc)) bundle$apparent_auc else NA_real_,
  train_oof_auc = if (!is.null(bundle$oof_auc)) bundle$oof_auc else NA_real_
)

fwrite(pred_ext_dt, out_pred_ext)
fwrite(perf_ext_dt, out_perf_ext)

png(out_roc_ext, width = 1800, height = 1600, res = 220)
plot(
  roc_ext,
  main = sprintf(
    "External ROC %s %s (%s, AUC = %.3f)",
    CFG$test_gse_id,
    TEST_REGION_USE,
    "AD vs Control",
    auc_ext
  ),
  legacy.axes = TRUE
)
abline(a = 0, b = 1, lty = 2)
dev.off()

# ============================================================
# 5) Store in state
# ============================================================
state$external_prob <- prob_ext
state$external_pred_class <- pred_class
state$external_auc <- auc_ext
state$external_roc <- roc_ext
state$external_perf <- perf_ext_dt
state$external_pred_dt <- pred_ext_dt
state$external_X_test_used <- X_test
state$external_missing_model_genes <- missing_cols


log_msg(
  "DONE EXTERNAL VALIDATION. dataset=", CFG$test_gse_id,
  " region=", TEST_REGION_USE,
  " n=", nrow(X_test),
  " AUC=", round(auc_ext, 4),
  " ACC@0.5=", round(acc_ext, 4),
  " SEN@0.5=", round(sens_ext, 4),
  " SPE@0.5=", round(spec_ext, 4)
)
