# Central configuration (edit this file)

# -------------------------
# Input files
# -------------------------
CORE_GENES_CSV <- "data/processed/core_genes.csv"                 # must contain column: GeneSymbol
GSE33000_SERIES_MATRIX_GZ <- "data/raw/GSE33000_series_matrix.txt.gz"

# optional validation dataset (set to NULL to skip)
VALID_SERIES_MATRIX_GZ <- "data/raw/GSE44770_series_matrix.txt.gz"

# If you have a curated metadata file, set it here; otherwise the pipeline will
# attempt to infer condition labels from GEO metadata (often unreliable).
GSE33000_SAMPLE_META_CSV <- "data/processed/sample_metadata_gse33000.csv"   # columns: sample_id, condition
VALID_SAMPLE_META_CSV    <- NULL

# -------------------------
# Platform annotation
# -------------------------
GPL_ID <- "GPL4372"  # platform for GSE33000 series matrix

# -------------------------
# DEG / effect filter
# -------------------------
FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.2

# -------------------------
# LASSO
# -------------------------
SEED <- 20260301
CV_FOLDS <- 10
USE_LAMBDA <- "lambda.1se"  # "lambda.min" or "lambda.1se"

# -------------------------
# Output folders
# -------------------------
DIR_DATA_PROCESSED <- "data/processed"
DIR_OUT_TABLES <- "outputs/tables"
DIR_OUT_FIGURES <- "outputs/figures"
DIR_OUT_LOGS <- "outputs/logs"
