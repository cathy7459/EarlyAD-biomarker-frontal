# ============================
# Configuration (edit this)
# ============================
CFG <- list(
  # --- TRAIN (GSE33000) ---
  train_series_matrix_gz = "data/raw/GSE33000_series_matrix.txt.gz",
  train_platform_id      = "GPL4372",
  core_gene_csv          = "outputs/WGCNA_significant_genes.csv",

  # --- TEST (GSE122063) ---
  test_gse_id            = "GSE122063",
  test_diag_col          = "patient diagnosis:ch1",

  # --- MODEL ---
  alpha                  = 1.0,     # 1=LASSO, 0=Ridge
  nfolds                 = 10,
  lambda_choice          = "lambda.1se",  # "lambda.min" or "lambda.1se"
  model_bundle_rds = "outputs/trainGSE33000_testGSE122063_model_bundle.rds",
  
  # --- Output ---
  out_dir                = "outputs",

)


# ------------------------------------------------------------
# 0) Config
# ------------------------------------------------------------
if (!exists("CFG")) stop("CFG not found. source('config/config.R') first.")

if (is.null(CFG$test_gse_id) || CFG$test_gse_id == "") {
  stop("CFG$test_gse_id is missing.")
}

if (is.null(CFG$test_diag_col) || CFG$test_diag_col == "") {
  stop("CFG$test_diag_col is missing.")
}

# brain region column in GSE122063 pData
if (!exists("TEST_REGION_COL")) {
  TEST_REGION_COL <- "brain region:ch1"
}

out_dir_raw <- file.path(CFG$out_dir, "split_raw")
dir.create(out_dir_raw, recursive = TRUE, showWarnings = FALSE)
