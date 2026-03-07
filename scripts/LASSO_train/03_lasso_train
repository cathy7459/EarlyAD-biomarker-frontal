
# scripts/03_lasso_train.R
source(file.path("scripts", "00_setup.R"))

log_msg("RUNNING FILE: scripts/03_lasso_train_validate.R")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

# ------------------------------------------------------------
# 0) Paths
# ------------------------------------------------------------
expr_path_candidates <- c(
  file.path("data", "processed", "gse33000_expr_filtered_postproj_gene_x_sample.csv"),
  file.path("data", "processed", "gse33000_expr_filtered_gene_x_sample.csv")
)

expr_path <- expr_path_candidates[file.exists(expr_path_candidates)][1]

if (is.na(expr_path) || !nzchar(expr_path)) {
  stop(
    "Missing filtered expression matrix.\nChecked:\n",
    paste(expr_path_candidates, collapse = "\n"),
    "\nRun scripts/02_deg_effect_filter.R first."
  )
}

meta_path <- file.path("data", "processed", "sample_metadata_gse33000.csv")
if (!file.exists(meta_path)) {
  stop("Missing metadata: ", meta_path)
}

out_sel <- file.path("outputs", "tables", "lasso_selected_genes.csv")
out_coef <- file.path("outputs", "tables", "lasso_coefficients.csv")
out_pred <- file.path("outputs", "tables", "lasso_train_predictions.csv")
out_roc  <- file.path("outputs", "figures", "roc_train.png")

dir.create(dirname(out_sel), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_roc), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1) Load expression
# ------------------------------------------------------------
log_msg("[1/4] Loading filtered expression: ", expr_path)

expr_dt <- fread(expr_path)
stopifnot("GeneSymbol" %in% names(expr_dt))

expr <- as.data.frame(expr_dt[, -1, with = FALSE])
rownames(expr) <- expr_dt$GeneSymbol

# genes x samples -> ensure numeric
expr[] <- lapply(expr, as.numeric)

# ------------------------------------------------------------
# 2) Load metadata and align samples
# ------------------------------------------------------------
log_msg("[2/4] Loading metadata: ", meta_path)

meta <- as.data.frame(fread(meta_path))

required_meta_cols <- c("sample_id", "condition")
if (!all(required_meta_cols %in% names(meta))) {
  stop("Metadata must contain columns: ", paste(required_meta_cols, collapse = ", "))
}

meta$sample_id <- trimws(as.character(meta$sample_id))
meta$condition <- trimws(as.character(meta$condition))

# keep only AD / Control
meta <- meta[meta$condition %in% c("AD", "Control"), , drop = FALSE]
meta <- meta[!duplicated(meta$sample_id), , drop = FALSE]

common_samples <- intersect(colnames(expr), meta$sample_id)
if (length(common_samples) < 20) {
  stop("Too few matched samples between expression and metadata: ", length(common_samples))
}

expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), , drop = FALSE]

stopifnot(identical(colnames(expr), meta$sample_id))

log_msg("Matched samples: ", ncol(expr))
log_msg("AD samples: ", sum(meta$condition == "AD"))
log_msg("Control samples: ", sum(meta$condition == "Control"))

# ------------------------------------------------------------
# 3) Build glmnet input
# ------------------------------------------------------------
# glmnet wants samples x genes
x <- t(as.matrix(expr))
storage.mode(x) <- "double"

if (anyNA(x)) {
  
  log_msg("NA detected in x -> performing median imputation")
  
  for (j in seq_len(ncol(x))) {
    idx <- is.na(x[,j])
    if (any(idx)) {
      x[idx,j] <- median(x[,j], na.rm=TRUE)
    }
  }
  
  log_msg("Remaining NA: ", sum(is.na(x)))
}

# outcome: Control = 0, AD = 1
y <- factor(meta$condition, levels = c("Control", "AD"))
if (anyNA(y)) stop("Outcome y contains NA.")
y_bin <- ifelse(y == "AD", 1, 0)

# remove zero-variance features
feature_sd <- apply(x, 2, sd, na.rm = TRUE)
keep_feat <- is.finite(feature_sd) & !is.na(feature_sd) & feature_sd > 0

if (sum(keep_feat) == 0) {
  stop("No non-constant features available for LASSO.")
}

x <- x[, keep_feat, drop = FALSE]

# final checks
if (anyNA(x)) {
  na_idx <- which(is.na(x), arr.ind = TRUE)
  print(head(na_idx, 20))
  stop("x has missing values.")
}
if (any(!is.finite(x))) {
  bad_idx <- which(!is.finite(x), arr.ind = TRUE)
  print(head(bad_idx, 20))
  stop("x has non-finite values (Inf / -Inf / NaN).")
}

log_msg("[3/4] LASSO input ready: samples=", nrow(x), " genes=", ncol(x))

# Save train statistics for external validation
feature_medians <- apply(x, 2, median, na.rm = TRUE)
feature_means <- colMeans(x, na.rm = TRUE)
feature_sds <- apply(x, 2, sd, na.rm = TRUE)

# safety
feature_medians[!is.finite(feature_medians)] <- 0
feature_means[!is.finite(feature_means)] <- 0
feature_sds[!is.finite(feature_sds) | feature_sds == 0] <- 1


# ------------------------------------------------------------
# 4) Train cv.glmnet
# ------------------------------------------------------------
if (!exists("SEED")) SEED <- 1234
if (!exists("CV_FOLDS")) CV_FOLDS <- 10
if (!exists("USE_LAMBDA")) USE_LAMBDA <- "lambda.1se"  # or "lambda.min"

set.seed(SEED)

log_msg("[4/4] Training LASSO (cv.glmnet) ...")
cvfit <- cv.glmnet(
  x = x,
  y = y_bin,
  family = "binomial",
  alpha = 1,
  nfolds = CV_FOLDS,
  standardize = TRUE,
  type.measure = "auc"
)

lam <- switch(
  USE_LAMBDA,
  "lambda.min" = cvfit$lambda.min,
  "lambda.1se" = cvfit$lambda.1se,
  stop("USE_LAMBDA must be 'lambda.min' or 'lambda.1se'")
)

log_msg("Selected lambda mode: ", USE_LAMBDA)
log_msg("Selected lambda value: ", signif(lam, 4))

# ------------------------------------------------------------
# 5) Extract selected genes / coefficients
# ------------------------------------------------------------
coef_mat <- as.matrix(coef(cvfit, s = lam))
coef_dt <- data.table(
  Feature = rownames(coef_mat),
  Coefficient = as.numeric(coef_mat[, 1])
)

# remove intercept
coef_dt_noint <- coef_dt[Feature != "(Intercept)"]

sel_dt <- coef_dt_noint[Coefficient != 0]

sel_dt[, abs_coef := abs(Coefficient)]
setorder(sel_dt, -abs_coef)
sel_dt[, abs_coef := NULL]

fwrite(sel_dt, out_sel)
fwrite(coef_dt, out_coef)

log_msg("Saved selected genes: ", out_sel)
log_msg("Saved full coefficient table: ", out_coef)
log_msg("Selected nonzero genes: ", nrow(sel_dt))

# ------------------------------------------------------------
# 6) Train predictions / ROC
# ------------------------------------------------------------
prob_train <- as.numeric(predict(cvfit, newx = x, s = lam, type = "response"))

pred_dt <- data.table(
  sample_id = rownames(x),
  condition = as.character(y),
  y_bin = y_bin,
  prob_AD = prob_train
)
fwrite(pred_dt, out_pred)
log_msg("Saved train predictions: ", out_pred)

roc_obj <- pROC::roc(response = y_bin, predictor = prob_train, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

png(out_roc, width = 1800, height = 1600, res = 220)
plot(
  roc_obj,
  main = sprintf("Training ROC (AUC = %.3f)", auc_val),
  legacy.axes = TRUE
)
abline(a = 0, b = 1, lty = 2)
dev.off()

log_msg("Saved ROC figure: ", out_roc)
log_msg("Training AUC: ", sprintf("%.3f", auc_val))

log_msg("DONE. samples=", nrow(x), " genes=", ncol(x), " selected=", nrow(sel_dt))

# ------------------------------------------------------------
# 7) Save model bundle for external validation
# ------------------------------------------------------------
bundle_prefix <- if (exists("CFG") && !is.null(CFG$prefix)) CFG$prefix else "ad_lasso"
out_bundle <- file.path("outputs", paste0(bundle_prefix, "_model_bundle.rds"))

bundle <- list(
  created_at = Sys.time(),
  script = "03_lasso_train_validate.R",
  
  cvfit = cvfit,
  lambda_use = lam,
  lambda_mode = USE_LAMBDA,
  
  genes_used = colnames(x),
  selected = sel_dt,
  coef_table = coef_dt,
  
  feature_medians = feature_medians,
  feature_means = feature_means,
  feature_sds = feature_sds,
  
  n_samples = nrow(x),
  n_genes = ncol(x),
  sample_ids = rownames(x),
  y_train = y_bin,
  y_train_chr = as.character(y),
  
  seed = SEED,
  cv_folds = CV_FOLDS
)

saveRDS(bundle, out_bundle)
log_msg("Saved model bundle: ", out_bundle)
