log_msg("RUNNING FILE: scripts/05_external_validate_gse122063.R")

if (!exists("state")) state <- new.env(parent = emptyenv())
if (!exists("CFG")) stop("CFG not found. source('config/config.R') first.")

# ------------------------------------------------------------
# 0) Check required objects from Step 04
# ------------------------------------------------------------
required_state_objects <- c("X_test", "y_test", "y_test_chr", "test_region_used")
missing_state_objects <- required_state_objects[!vapply(required_state_objects, exists, logical(1), envir = state)]

if (length(missing_state_objects) > 0) {
  stop(
    "Missing required objects from Step 04: ",
    paste(missing_state_objects, collapse = ", "),
    "\nRun scripts/04_prepare_external_test_set.R first."
  )
}

if (is.null(CFG$model_bundle_rds) || !nzchar(CFG$model_bundle_rds) || !file.exists(CFG$model_bundle_rds)) {
  stop("Model bundle not found: ", CFG$model_bundle_rds)
}

out_tables  <- file.path(CFG$out_dir, "tables")
out_figures <- file.path(CFG$out_dir, "figures")
.dir_create(out_tables)
.dir_create(out_figures)

# ------------------------------------------------------------
# 1) Recover model bundle and validation inputs
# ------------------------------------------------------------
bundle <- readRDS(CFG$model_bundle_rds)

if (is.null(bundle$cvfit)) {
  stop("cvfit missing in model bundle.")
}
if (is.null(bundle$cvfit$glmnet.fit)) {
  stop("glmnet.fit missing in bundle$cvfit.")
}
if (is.null(bundle$lambda_use) || !is.finite(bundle$lambda_use)) {
  stop("lambda_use missing or invalid in model bundle.")
}
if (is.null(bundle$genes_used) || length(bundle$genes_used) < 2) {
  stop("genes_used missing or too short in model bundle.")
}

TEST_REGION_USE <- state$test_region_used
region_tag <- paste0(TEST_REGION_USE, "_excluding_vascular")

bundle_prefix <- if (!is.null(CFG$prefix) && nzchar(CFG$prefix)) {
  CFG$prefix
} else {
  "ad_enet"
}

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

out_coef_path <- file.path(
  out_figures,
  paste0(bundle_prefix, "_model_", region_tag, "_coefficient_path.png")
)

model_fit <- bundle$cvfit$glmnet.fit
lambda_use <- bundle$lambda_use
genes_used <- as.character(bundle$genes_used)

X_test <- state$X_test

if (!all(genes_used %in% colnames(X_test))) {
  missing_model_genes <- setdiff(genes_used, colnames(X_test))
  stop(
    "Some model genes are missing in state$X_test: ",
    paste(head(missing_model_genes, 20), collapse = ", "),
    if (length(missing_model_genes) > 20) " ..." else ""
  )
}

X_test <- X_test[, genes_used, drop = FALSE]

# ------------------------------------------------------------
# 2) External prediction
# ------------------------------------------------------------
prob_ext <- as.numeric(
  predict(
    model_fit,
    newx = X_test,
    s = lambda_use,
    type = "response"
  )
)

roc_obj <- pROC::roc(
  response = state$y_test,
  predictor = prob_ext,
  quiet = TRUE
)

auc_ext <- as.numeric(pROC::auc(roc_obj))

pred_class <- ifelse(prob_ext >= 0.5, "AD", "Control")
pred_bin <- ifelse(pred_class == "AD", 1, 0)

acc_ext <- mean(pred_bin == state$y_test)
sens_ext <- if (sum(state$y_test == 1) > 0) {
  mean(pred_bin[state$y_test == 1] == 1)
} else {
  NA_real_
}
spec_ext <- if (sum(state$y_test == 0) > 0) {
  mean(pred_bin[state$y_test == 0] == 0)
} else {
  NA_real_
}

# ------------------------------------------------------------
# 3) Save prediction and performance tables
# ------------------------------------------------------------
pred_dt <- data.table::data.table(
  sample_id = rownames(X_test),
  condition = state$y_test_chr,
  y_bin = state$y_test,
  prob_AD = prob_ext,
  pred_class_0.5 = pred_class,
  pred_bin_0.5 = pred_bin
)

perf_dt <- data.table::data.table(
  dataset = CFG$test_gse_id,
  region = TEST_REGION_USE,
  auc = auc_ext,
  accuracy_0.5 = acc_ext,
  sensitivity_0.5 = sens_ext,
  specificity_0.5 = spec_ext,
  n_samples = nrow(X_test),
  n_AD = sum(state$y_test == 1),
  n_Control = sum(state$y_test == 0),
  alpha_use = if (!is.null(bundle$alpha_use)) bundle$alpha_use else NA_real_,
  lambda_use = lambda_use,
  matched_genes = if (!is.null(state$test_matched_genes)) length(state$test_matched_genes) else NA_integer_,
  missing_genes_imputed = if (!is.null(state$test_missing_genes)) length(state$test_missing_genes) else NA_integer_
)

data.table::fwrite(pred_dt, out_pred_ext)
data.table::fwrite(perf_dt, out_perf_ext)

# ------------------------------------------------------------
# 4) Save ROC figure
# ------------------------------------------------------------
png(out_roc_ext, width = 1800, height = 1600, res = 220)
plot(
  roc_obj,
  main = sprintf(
    "External ROC: %s (%s) | AUC = %.3f",
    CFG$test_gse_id,
    TEST_REGION_USE,
    auc_ext
  ),
  legacy.axes = TRUE
)
abline(a = 0, b = 1, lty = 2)
dev.off()

# ------------------------------------------------------------
# 5) Save coefficient path plot
# ------------------------------------------------------------
png(out_coef_path, width = 2000, height = 1700, res = 220)
plot(
  model_fit,
  xvar = "lambda",
  label = FALSE,
  main = sprintf(
    "Coefficient Path (%s model | alpha = %s)\n",
    TEST_REGION_USE,
    if (!is.null(bundle$alpha_use)) format(round(bundle$alpha_use, 3), nsmall = 3) else "NA"
  )
)
abline(v = log(lambda_use), lty = 2)
mtext(
  side = 3,
  line = 0.5,
  adj = 1,
  text = paste0("lambda_use = ", signif(lambda_use, 4)),
  cex = 0.9
)
dev.off()

# ------------------------------------------------------------
# 6) Store results in state
# ------------------------------------------------------------
state$external_prob <- prob_ext
state$external_pred_class_0.5 <- pred_class
state$external_pred_bin_0.5 <- pred_bin
state$external_roc <- roc_obj
state$external_auc <- auc_ext
state$external_perf <- perf_dt
state$external_coef_path_file <- out_coef_path

log_msg("[SUMMARY] External validation dataset: ", CFG$test_gse_id)
log_msg("[SUMMARY] Region: ", TEST_REGION_USE)
log_msg("[SUMMARY] Samples: ", nrow(X_test))
log_msg("[SUMMARY] AUC: ", round(auc_ext, 4))
log_msg("[SUMMARY] Accuracy (0.5): ", round(acc_ext, 4))
log_msg("[SUMMARY] Sensitivity (0.5): ", round(sens_ext, 4))
log_msg("[SUMMARY] Specificity (0.5): ", round(spec_ext, 4))
log_msg("[SUMMARY] ROC saved to: ", out_roc_ext)
log_msg("[SUMMARY] Coefficient path saved to: ", out_coef_path)

log_msg("DONE external validation. AUC = ", round(auc_ext, 4))

