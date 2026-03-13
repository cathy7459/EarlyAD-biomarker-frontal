# config/enet.R
# Central configuration for the Early AD Elastic Net pipeline.
# R version target: 4.5.1

CFG <- list(
  prefix = "ad_enet",

  # --- TRAIN (GSE33000) ---
  train_series_matrix_gz = "data/raw/GSE33000_series_matrix.txt.gz",
  train_platform_id      = "GPL4372",
  core_gene_csv          = "outputs/WGCNA_significant_genes.csv",

  # --- TEST (GSE122063) ---
  test_gse_id            = "GSE122063",
  test_diag_col          = "patient diagnosis:ch1",
  test_region            = "frontal",

  # --- MODEL / BUNDLE ---
  model_bundle_rds       = "outputs/ad_enet_model_bundle.rds",

  # --- Output ---
  out_dir                = "outputs"
)

# ============================================================
# Hyperparameters
# ============================================================
SEED <- 1234

# moderately relaxed, but still conservative
UNIV_TOP_N <- 80
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

# for DEG filter stage
LOGFC_CUTOFF <- 0.2
