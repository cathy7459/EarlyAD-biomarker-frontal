# R/utils_model.R
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

# ------------------------------------------------------------
# Correlation pruning
# Removes highly correlated features while keeping one feature
# from each correlated cluster.
# Input x must be a numeric matrix/data.frame with samples in rows
# and features in columns.
# ------------------------------------------------------------
prune_correlated_features <- function(x, cor_cutoff = 0.90) {
  # Coerce to numeric matrix
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if (is.null(colnames(x))) {
    stop("x must have column names for correlation pruning.")
  }
  
  if (ncol(x) <= 1) {
    return(colnames(x))
  }
  
  # Drop zero-variance or all-NA features first
  sds <- apply(x, 2, function(v) sd(v, na.rm = TRUE))
  keep <- !(is.na(sds) | sds == 0)
  
  if (!any(keep)) {
    stop("No usable features remained after removing zero-variance/all-NA columns.")
  }
  
  x_use <- x[, keep, drop = FALSE]
  
  if (ncol(x_use) <= 1) {
    return(colnames(x_use))
  }
  
  # Pairwise-complete correlation for robustness to residual NAs
  cor_mat <- suppressWarnings(stats::cor(
    x_use,
    use = "pairwise.complete.obs",
    method = "pearson"
  ))
  
  # Replace NA correlations with 0 so they are not spuriously removed
  cor_mat[is.na(cor_mat)] <- 0
  
  # Work on absolute correlation and ignore diagonal
  abs_cor <- abs(cor_mat)
  diag(abs_cor) <- 0
  
  feature_names <- colnames(x_use)
  remove_idx <- rep(FALSE, length(feature_names))
  
  # Greedy pruning:
  # for each highly correlated pair, remove the feature with
  # the larger mean absolute correlation to all other features
  repeat {
    if (ncol(abs_cor) <= 1) break
    
    max_cor <- max(abs_cor, na.rm = TRUE)
    if (!is.finite(max_cor) || max_cor < cor_cutoff) break
    
    hit <- which(abs_cor == max_cor, arr.ind = TRUE)[1, ]
    i <- hit[1]
    j <- hit[2]
    
    mean_i <- mean(abs_cor[i, ], na.rm = TRUE)
    mean_j <- mean(abs_cor[j, ], na.rm = TRUE)
    
    drop_pos <- if (mean_i >= mean_j) i else j
    drop_name <- colnames(abs_cor)[drop_pos]
    
    remove_idx[match(drop_name, feature_names)] <- TRUE
    
    keep_now <- !colnames(abs_cor) %in% drop_name
    abs_cor <- abs_cor[keep_now, keep_now, drop = FALSE]
  }
  
  kept_features <- feature_names[!remove_idx]
  
  if (length(kept_features) == 0) {
    stop("Correlation pruning removed all features.")
  }
  
  kept_features
}

make_stratified_foldid <- function(y_bin, k = 5, seed = 1) {
  set.seed(seed)
  foldid <- integer(length(y_bin))

  idx0 <- which(y_bin == 0)
  idx1 <- which(y_bin == 1)

  foldid[idx0] <- sample(rep(seq_len(k), length.out = length(idx0)))
  foldid[idx1] <- sample(rep(seq_len(k), length.out = length(idx1)))

  foldid
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

  list(ranking = univ_dt, selected_features = selected)
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

save_roc_plot <- function(roc_obj, file, title) {
  png(file, width = 1800, height = 1600, res = 220)
  plot(roc_obj, main = title, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2)
  dev.off()
}
