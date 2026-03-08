# scripts/03b_enet_train_validate.R

log_msg("RUNNING FILE: scripts/03b_enet_train_validate.R")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

if (!exists("state")) state <- new.env(parent = emptyenv())
if (!exists("CFG")) stop("CFG not found. source('config/config.R') first.")

# ============================================================
# hyperparameter
# ============================================================
SEED <- 1234

# moderately relaxed, but still conservative
UNIV_TOP_N <- 60
UNIV_METHOD <- "wilcox"
USE_FDR_GATE <- TRUE
FDR_CUTOFF <- 0.30
MIN_FDR_PASS <- 8

CV_FOLDS <- 5
CV_REPEATS <- 30

# allow mild ENet relaxation from pure LASSO
ALPHA_GRID <- c(0.40, 0.55, 0.70, 0.85, 1.00)
ALPHA_TOL <- 0.02
ALPHA_COMPLEXITY_PENALTY <- 0

# lambda is FIXED to 1se
LAMBDA_MODE <- "lambda.1se"

# keep final model reasonably sparse, but not too tiny
MAX_FINAL_GENES <- 25
MIN_FINAL_GENES <- 4

log_msg("[HP] lambda rule fixed to lambda.1se")
log_msg("[HP] UNIV_TOP_N=", UNIV_TOP_N,
        " | FDR gate=", USE_FDR_GATE,
        " | FDR_CUTOFF=", FDR_CUTOFF,
        " | CV_FOLDS=", CV_FOLDS,
        " | CV_REPEATS=", CV_REPEATS)

# ============================================================
# Paths
# ============================================================
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

out_tables <- file.path(CFG$out_dir, "tables")
out_figures <- file.path(CFG$out_dir, "figures")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)

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

if (is.null(CFG$model_bundle_rds) || !nzchar(CFG$model_bundle_rds)) {
  stop("CFG$model_bundle_rds is missing.")
}
out_bundle <- CFG$model_bundle_rds
dir.create(dirname(out_bundle), recursive = TRUE, showWarnings = FALSE)
# ============================================================
# Helper functions
# ============================================================
make_stratified_foldid <- function(y_bin, k = 5, seed = 1) {
  set.seed(seed)
  foldid <- integer(length(y_bin))
  
  idx0 <- which(y_bin == 0)
  idx1 <- which(y_bin == 1)
  
  foldid[idx0] <- sample(rep(seq_len(k), length.out = length(idx0)))
  foldid[idx1] <- sample(rep(seq_len(k), length.out = length(idx1)))
  
  foldid
}

extract_coef_table <- function(glmnet_fit, lambda) {
  coef_mat <- as.matrix(coef(glmnet_fit, s = lambda))
  data.table(
    Feature = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat[, 1])
  )
}

extract_selected_genes <- function(coef_dt) {
  sel_dt <- coef_dt[Feature != "(Intercept)" & Coefficient != 0]
  if (nrow(sel_dt) > 0) {
    sel_dt[, abs_coef := abs(Coefficient)]
    setorder(sel_dt, -abs_coef)
    sel_dt[, abs_coef := NULL]
  }
  sel_dt
}

run_univariate_filter <- function(x, y_bin,
                                  method = "wilcox",
                                  use_fdr_gate = TRUE,
                                  fdr_cutoff = 0.25,
                                  min_fdr_pass = 10,
                                  top_n = 50) {
  pvals <- rep(NA_real_, ncol(x))
  effs  <- rep(NA_real_, ncol(x))
  
  names(pvals) <- colnames(x)
  names(effs)  <- colnames(x)
  
  for (j in seq_len(ncol(x))) {
    xx <- x[, j]
    x_ad   <- xx[y_bin == 1]
    x_ctrl <- xx[y_bin == 0]
    
    pvals[j] <- switch(
      method,
      "wilcox" = tryCatch(wilcox.test(x_ad, x_ctrl, exact = FALSE)$p.value, error = function(e) 1),
      "ttest"  = tryCatch(t.test(x_ad, x_ctrl)$p.value, error = function(e) 1),
      stop("UNIV_METHOD must be 'wilcox' or 'ttest'")
    )
    
    effs[j] <- abs(median(x_ad, na.rm = TRUE) - median(x_ctrl, na.rm = TRUE))
  }
  
  pvals[!is.finite(pvals) | is.na(pvals)] <- 1
  effs[!is.finite(effs) | is.na(effs)] <- 0
  
  padj <- p.adjust(pvals, method = "BH")
  
  univ_dt <- data.table(
    Feature = colnames(x),
    p_value = pvals,
    p_adj = padj,
    effect_size = effs
  )
  
  if (use_fdr_gate) {
    gated_dt <- univ_dt[p_adj <= fdr_cutoff]
    if (nrow(gated_dt) >= min_fdr_pass) {
      setorder(gated_dt, p_adj, -effect_size)
      selected <- gated_dt$Feature[seq_len(min(top_n, nrow(gated_dt)))]
    } else {
      setorder(univ_dt, p_value, -effect_size)
      selected <- univ_dt$Feature[seq_len(min(top_n, nrow(univ_dt)))]
    }
  } else {
    setorder(univ_dt, p_value, -effect_size)
    selected <- univ_dt$Feature[seq_len(min(top_n, nrow(univ_dt)))]
  }
  
  list(
    ranking = univ_dt,
    selected_features = selected
  )
}

fit_cv_glmnet <- function(x, y_bin, alpha, foldid) {
  cv.glmnet(
    x = x,
    y = y_bin,
    family = "binomial",
    alpha = alpha,
    foldid = foldid,
    standardize = TRUE,
    type.measure = "auc"
  )
}

get_lambda_1se <- function(cvfit) {
  cvfit$lambda.1se
}

get_oof_predictions <- function(x, y_bin, alpha, foldid, lambda_value) {
  n <- nrow(x)
  pred <- rep(NA_real_, n)
  
  for (k in sort(unique(foldid))) {
    idx_te <- which(foldid == k)
    idx_tr <- which(foldid != k)
    
    fit_k <- glmnet(
      x = x[idx_tr, , drop = FALSE],
      y = y_bin[idx_tr],
      family = "binomial",
      alpha = alpha,
      lambda = lambda_value,
      standardize = TRUE
    )
    
    pred[idx_te] <- as.numeric(
      predict(fit_k, newx = x[idx_te, , drop = FALSE], s = lambda_value, type = "response")
    )
  }
  
  pred
}

# ============================================================
# 1) Load expression matrix
# ============================================================
log_msg("[1/7] Loading expression matrix: ", expr_path)

expr_dt <- fread(expr_path)
stopifnot("GeneSymbol" %in% names(expr_dt))

expr <- as.data.frame(expr_dt[, -1, with = FALSE])
rownames(expr) <- expr_dt$GeneSymbol
expr[] <- lapply(expr, as.numeric)

# optional core-gene restriction from prior step
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
  } else {
    warning("Too few overlapping core genes found; skipping core-gene restriction.")
  }
}

# ============================================================
# 2) Load metadata and align
# ============================================================
log_msg("[2/7] Loading metadata: ", meta_path)

meta <- as.data.frame(fread(meta_path))

required_meta_cols <- c("sample_id", "condition")
if (!all(required_meta_cols %in% names(meta))) {
  stop("Metadata must contain columns: ", paste(required_meta_cols, collapse = ", "))
}

meta$sample_id <- trimws(as.character(meta$sample_id))
meta$condition <- trimws(as.character(meta$condition))
meta <- meta[meta$condition %in% c("AD", "Control"), , drop = FALSE]
meta <- meta[!duplicated(meta$sample_id), , drop = FALSE]

common_samples <- intersect(colnames(expr), meta$sample_id)
if (length(common_samples) < 20) {
  stop("Too few matched samples between expression and metadata: ", length(common_samples))
}

expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
stopifnot(identical(colnames(expr), meta$sample_id))

log_msg("[2/7] Matched samples: ", ncol(expr))
log_msg("[2/7] AD samples: ", sum(meta$condition == "AD"))
log_msg("[2/7] Control samples: ", sum(meta$condition == "Control"))

# ============================================================
# 3) Build model input
# ============================================================
log_msg("[3/7] Building model matrix")

x <- t(as.matrix(expr))
storage.mode(x) <- "double"

if (anyNA(x)) {
  log_msg("[3/7] NA detected -> median imputation")
  for (j in seq_len(ncol(x))) {
    idx <- is.na(x[, j])
    if (any(idx)) {
      x[idx, j] <- median(x[, j], na.rm = TRUE)
    }
  }
}

y <- factor(meta$condition, levels = c("Control", "AD"))
if (anyNA(y)) stop("Outcome contains NA values.")
y_bin <- ifelse(y == "AD", 1, 0)

feature_sd <- apply(x, 2, sd, na.rm = TRUE)
keep_feat <- is.finite(feature_sd) & !is.na(feature_sd) & feature_sd > 0
x <- x[, keep_feat, drop = FALSE]

if (ncol(x) < 2) stop("Too few non-constant features after filtering.")

feature_medians_full <- apply(x, 2, median, na.rm = TRUE)
feature_means_full   <- colMeans(x, na.rm = TRUE)
feature_sds_full     <- apply(x, 2, sd, na.rm = TRUE)

feature_medians_full[!is.finite(feature_medians_full)] <- 0
feature_means_full[!is.finite(feature_means_full)] <- 0
feature_sds_full[!is.finite(feature_sds_full) | feature_sds_full == 0] <- 1

log_msg("[3/7] Input ready: samples=", nrow(x), ", genes=", ncol(x))

# ============================================================
# 4) Univariate prefilter
# ============================================================
log_msg("[4/7] Running univariate prefilter")

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

if (length(prefilter_genes) < 5) {
  stop("Too few genes selected by univariate prefilter.")
}

x_filt <- x[, prefilter_genes, drop = FALSE]
fwrite(univ_dt, out_univ)

log_msg("[4/7] Prefiltered genes: ", ncol(x_filt))

# ============================================================
# 5) Repeated CV for alpha selection
# ============================================================
log_msg("[5/7] Repeated stratified CV for alpha selection")

set.seed(SEED)
perf_list <- vector("list", length(ALPHA_GRID) * CV_REPEATS)
kk <- 1

for (a in ALPHA_GRID) {
  for (r in seq_len(CV_REPEATS)) {
    foldid_r <- make_stratified_foldid(
      y_bin = y_bin,
      k = CV_FOLDS,
      seed = SEED + 1000 * r + round(a * 1000)
    )
    
    cv_tmp <- fit_cv_glmnet(
      x = x_filt,
      y_bin = y_bin,
      alpha = a,
      foldid = foldid_r
    )
    
    lam_r <- get_lambda_1se(cv_tmp)
    coef_r <- extract_coef_table(cv_tmp$glmnet.fit, lam_r)
    sel_r <- extract_selected_genes(coef_r)
    n_nonzero_r <- nrow(sel_r)
    
    # use cv AUC estimate, not apparent AUC
    idx_1se <- which.min(abs(cv_tmp$lambda - lam_r))
    cv_auc_1se <- cv_tmp$cvm[idx_1se]
    cv_auc_1se_sd <- cv_tmp$cvsd[idx_1se]
    
    perf_list[[kk]] <- data.table(
      alpha = a,
      repeat_id = r,
      lambda_1se = lam_r,
      cv_auc_1se = cv_auc_1se,
      cv_auc_1se_sd = cv_auc_1se_sd,
      nonzero_genes = n_nonzero_r
    )
    kk <- kk + 1
  }
}

perf_dt <- rbindlist(perf_list)
fwrite(perf_dt, out_alpha_raw)

alpha_summary <- perf_dt[, .(
  mean_cv_auc_1se = mean(cv_auc_1se, na.rm = TRUE),
  sd_cv_auc_1se   = sd(cv_auc_1se, na.rm = TRUE),
  mean_nonzero    = mean(nonzero_genes, na.rm = TRUE),
  median_nonzero  = median(nonzero_genes, na.rm = TRUE),
  min_nonzero     = min(nonzero_genes, na.rm = TRUE),
  max_nonzero     = max(nonzero_genes, na.rm = TRUE)
), by = alpha]

alpha_summary[, conservative_score := mean_cv_auc_1se - sd_cv_auc_1se]
best_conservative_score <- max(alpha_summary$conservative_score, na.rm = TRUE)

# keep alpha candidates whose performance is practically equivalent
candidate_alpha <- alpha_summary[
  conservative_score >= (best_conservative_score - ALPHA_TOL)
]

# apply only a very light complexity penalty
candidate_alpha[, penalized_score := conservative_score - ALPHA_COMPLEXITY_PENALTY * mean_nonzero]

# among near-tied models, prefer smaller alpha (more stable ENet) over pure LASSO
setorder(candidate_alpha, -penalized_score, alpha, mean_nonzero)

best_alpha <- candidate_alpha$alpha[1]
fwrite(alpha_summary, out_alpha)

log_msg("[5/7] Alpha summary:")
print(alpha_summary)

log_msg("[5/7] Alpha candidates within tolerance:")
print(candidate_alpha)

log_msg("[5/7] Selected alpha: ", best_alpha)
# ============================================================
# 6) Final fit + OOF evaluation
# ============================================================
log_msg("[6/7] Final model fit with selected alpha")

final_foldid <- make_stratified_foldid(
  y_bin = y_bin,
  k = CV_FOLDS,
  seed = SEED + 99999
)

cvfit <- fit_cv_glmnet(
  x = x_filt,
  y_bin = y_bin,
  alpha = best_alpha,
  foldid = final_foldid
)

lam <- get_lambda_1se(cvfit)
coef_dt <- extract_coef_table(cvfit$glmnet.fit, lam)
sel_dt <- extract_selected_genes(coef_dt)

if (nrow(sel_dt) > MAX_FINAL_GENES) {
  warning(
    "Selected genes (", nrow(sel_dt), ") exceed MAX_FINAL_GENES (", MAX_FINAL_GENES,
    "). This is still acceptable because lambda.1se is fixed and should not be tightened further."
  )
}
if (nrow(sel_dt) < MIN_FINAL_GENES) {
  warning("Too few genes selected; model may still be overly strict.")
}

fwrite(sel_dt, out_sel)
fwrite(coef_dt, out_coef)

# Apparent prediction (for reference only)
prob_app <- as.numeric(
  predict(cvfit$glmnet.fit, newx = x_filt, s = lam, type = "response")
)

pred_app_dt <- data.table(
  sample_id = rownames(x_filt),
  condition = as.character(y),
  y_bin = y_bin,
  prob_AD = prob_app
)
fwrite(pred_app_dt, out_pred_app)

roc_app <- pROC::roc(response = y_bin, predictor = prob_app, quiet = TRUE)
auc_app <- as.numeric(pROC::auc(roc_app))

png(out_roc_app, width = 1800, height = 1600, res = 220)
plot(
  roc_app,
  main = sprintf("Training Apparent ROC Elastic Net (AUC = %.3f)", auc_app),
  legacy.axes = TRUE
)
abline(a = 0, b = 1, lty = 2)
dev.off()

# OOF prediction (more realistic training estimate)
prob_oof <- get_oof_predictions(
  x = x_filt,
  y_bin = y_bin,
  alpha = best_alpha,
  foldid = final_foldid,
  lambda_value = lam
)

pred_oof_dt <- data.table(
  sample_id = rownames(x_filt),
  condition = as.character(y),
  y_bin = y_bin,
  prob_AD = prob_oof
)
fwrite(pred_oof_dt, out_pred_oof)

roc_oof <- pROC::roc(response = y_bin, predictor = prob_oof, quiet = TRUE)
auc_oof <- as.numeric(pROC::auc(roc_oof))

png(out_roc_oof, width = 1800, height = 1600, res = 220)
plot(
  roc_oof,
  main = sprintf("Training OOF ROC Elastic Net (AUC = %.3f)", auc_oof),
  legacy.axes = TRUE
)
abline(a = 0, b = 1, lty = 2)
dev.off()

# ============================================================
# 7) Save model bundle
# ============================================================
bundle <- list(
  created_at = Sys.time(),
  script = "03b_enet_train_validate.R",
  
  cvfit = cvfit,
  lambda_use = lam,
  lambda_mode = LAMBDA_MODE,
  alpha_use = best_alpha,
  
  genes_used = colnames(x_filt),        # exact test alignment space
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

# store in state for downstream scripts if needed
state$train_x_filt <- x_filt
state$train_y_bin <- y_bin
state$train_meta <- meta
state$enet_bundle <- bundle

log_msg(
  "DONE. samples=", nrow(x_filt),
  " genes=", ncol(x_filt),
  " selected=", nrow(sel_dt),
  " apparent_AUC=", round(auc_app, 4),
  " OOF_AUC=", round(auc_oof, 4),
  " bundle=", out_bundle
)
