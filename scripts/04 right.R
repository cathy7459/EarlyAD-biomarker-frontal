log_msg("=== PREPARE EXTERNAL TEST SET FROM PRE-SPLIT OBJECTS: ", CFG$test_gse_id, " ===")

if (!exists("state")) state <- new.env(parent = emptyenv())

if (is.null(CFG$model_bundle_rds) || !nzchar(CFG$model_bundle_rds) || !file.exists(CFG$model_bundle_rds)) {
  stop("Model bundle not found: ", CFG$model_bundle_rds)
}

.dir_create(file.path(CFG$out_dir, "tables"))

# ------------------------------------------------------------
# 0) Select region
# ------------------------------------------------------------
TEST_REGION_USE <- if (!is.null(CFG$test_region) && nzchar(CFG$test_region)) {
  tolower(trimws(CFG$test_region))
} else {
  "frontal"
}

if (!TEST_REGION_USE %in% c("frontal", "temporal")) {
  stop("CFG$test_region must be either 'frontal' or 'temporal'. Got: ", TEST_REGION_USE)
}

# ------------------------------------------------------------
# 1) Load model bundle and recover training-side feature stats
# ------------------------------------------------------------
bundle <- readRDS(CFG$model_bundle_rds)

if (is.null(bundle$genes_used)) {
  stop("Model bundle does not contain 'genes_used'.")
}
genes_used <- unique(as.character(bundle$genes_used))
genes_used <- genes_used[!is.na(genes_used) & nzchar(genes_used)]

if (length(genes_used) == 0) {
  stop("No valid genes found in bundle$genes_used.")
}

# Helper: recover named numeric vector from multiple possible bundle keys
.get_named_numeric_from_bundle <- function(bundle, candidates, required = TRUE, fill_value = NA_real_) {
  obj <- NULL
  hit_name <- NA_character_
  
  for (nm in candidates) {
    if (!is.null(bundle[[nm]])) {
      obj <- bundle[[nm]]
      hit_name <- nm
      break
    }
  }
  
  if (is.null(obj)) {
    if (required) {
      stop("None of the required fields were found in model bundle: ",
           paste(candidates, collapse = ", "))
    } else {
      out <- setNames(rep(fill_value, length(genes_used)), genes_used)
      return(out)
    }
  }
  
  obj <- unlist(obj, use.names = TRUE)
  
  if (is.null(names(obj))) {
    stop("Bundle field found (", hit_name, ") but it does not have gene names.")
  }
  
  obj <- as.numeric(obj)
  names(obj) <- names(unlist(bundle[[hit_name]], use.names = TRUE))
  
  out <- obj[genes_used]
  names(out) <- genes_used
  out
}

# Training medians for missing-value imputation
train_gene_medians <- .get_named_numeric_from_bundle(
  bundle,
  candidates = c("feature_medians", "train_feature_medians", "medians", "feature_median"),
  required = FALSE,
  fill_value = 0
)
train_gene_medians[is.na(train_gene_medians) | !is.finite(train_gene_medians)] <- 0

# Training means for scaling
train_gene_means <- .get_named_numeric_from_bundle(
  bundle,
  candidates = c("feature_means", "feature_center", "train_means", "center"),
  required = TRUE
)

# Training SDs for scaling
train_gene_sds <- .get_named_numeric_from_bundle(
  bundle,
  candidates = c("feature_sds", "feature_scale", "train_sds", "scale"),
  required = TRUE
)

train_gene_means[is.na(train_gene_means) | !is.finite(train_gene_means)] <- 0
train_gene_sds[is.na(train_gene_sds) | !is.finite(train_gene_sds) | train_gene_sds <= 0] <- 1

# ------------------------------------------------------------
# 2) Load pre-split test objects
# ------------------------------------------------------------
if (TEST_REGION_USE == "frontal") {
  expr_gene_gxs <- state$gse122063_frontal_expr_gene_gxs
  meta_sub <- state$gse122063_frontal_meta
  expr_probe_sub <- state$gse122063_frontal_expr_probe
} else {
  expr_gene_gxs <- state$gse122063_temporal_expr_gene_gxs
  meta_sub <- state$gse122063_temporal_meta
  expr_probe_sub <- state$gse122063_temporal_expr_probe
}

if (is.null(expr_gene_gxs) || is.null(meta_sub) || is.null(expr_probe_sub)) {
  stop(
    "Pre-split test objects are missing for region '", TEST_REGION_USE, "'. ",
    "Run the GSE122063 split/collapse step first."
  )
}

if (!(CFG$test_diag_col %in% colnames(meta_sub))) {
  hit <- colnames(meta_sub)[tolower(colnames(meta_sub)) == tolower(CFG$test_diag_col)]
  if (length(hit) == 1) {
    CFG$test_diag_col <- hit
  } else {
    stop("Diagnosis column not found in test metadata: ", CFG$test_diag_col)
  }
}

# ------------------------------------------------------------
# 3) Enforce sample alignment between expression and metadata
# ------------------------------------------------------------
expr_samples <- colnames(expr_gene_gxs)
meta_samples <- rownames(meta_sub)

if (is.null(expr_samples) || length(expr_samples) == 0) {
  stop("Expression matrix has no sample column names.")
}
if (is.null(meta_samples) || length(meta_samples) == 0) {
  stop("Metadata has no row names for sample IDs.")
}

common_samples <- intersect(expr_samples, meta_samples)
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between expression matrix and metadata.")
}

expr_gene_gxs <- expr_gene_gxs[, common_samples, drop = FALSE]
expr_probe_sub <- expr_probe_sub[, common_samples, drop = FALSE]
meta_sub <- meta_sub[common_samples, , drop = FALSE]

# Safety check
if (!identical(colnames(expr_gene_gxs), rownames(meta_sub))) {
  stop("Sample alignment failed: expression columns and metadata rownames are not identical.")
}
if (!identical(colnames(expr_probe_sub), rownames(meta_sub))) {
  stop("Probe-level matrix and metadata are not aligned after subsetting.")
}

# ------------------------------------------------------------
# 4) Construct test labels (AD vs Control only)
# ------------------------------------------------------------
diag_sub_lc <- trimws(tolower(as.character(meta_sub[[CFG$test_diag_col]])))

valid_idx <- diag_sub_lc %in% c("alzheimer's disease", "control")
if (!all(valid_idx)) {
  log_msg("[WARN] Dropping non-AD/control samples at test stage: ", sum(!valid_idx))
  expr_gene_gxs <- expr_gene_gxs[, valid_idx, drop = FALSE]
  expr_probe_sub <- expr_probe_sub[, valid_idx, drop = FALSE]
  meta_sub <- meta_sub[valid_idx, , drop = FALSE]
  diag_sub_lc <- diag_sub_lc[valid_idx]
}

if (ncol(expr_gene_gxs) == 0) {
  stop("No valid AD/control samples remain in the external test set.")
}

state$test_samples <- colnames(expr_probe_sub)
state$y_test_chr <- ifelse(diag_sub_lc == "alzheimer's disease", "AD", "Control")
state$y_test <- ifelse(state$y_test_chr == "AD", 1, 0)

# ------------------------------------------------------------
# 4.5) Apply the same training-side projected feature space
# ------------------------------------------------------------
train_proj_path <- file.path("data", "processed", "gse33000_expr_projected_gene_x_sample.csv")
if (!file.exists(train_proj_path)) {
  stop("Training projected feature-space file not found: ", train_proj_path)
}

train_proj_dt <- data.table::fread(train_proj_path)
if (!("GeneSymbol" %in% names(train_proj_dt))) {
  stop("Projected training matrix must contain column 'GeneSymbol': ", train_proj_path)
}

projected_genes <- unique(trimws(as.character(train_proj_dt$GeneSymbol)))
projected_genes <- projected_genes[!is.na(projected_genes) & nzchar(projected_genes)]

test_gene_mat_full <- as.matrix(expr_gene_gxs)
storage.mode(test_gene_mat_full) <- "double"

if (is.null(rownames(test_gene_mat_full))) {
  stop("Gene-level test matrix has no rownames (gene symbols).")
}

matched_projected_genes <- intersect(projected_genes, rownames(test_gene_mat_full))

if (length(matched_projected_genes) == 0) {
  stop("No overlap between training projected genes and external test gene matrix.")
}

test_gene_mat <- matrix(
  NA_real_,
  nrow = length(projected_genes),
  ncol = ncol(test_gene_mat_full),
  dimnames = list(projected_genes, colnames(test_gene_mat_full))
)

test_gene_mat[matched_projected_genes, ] <- test_gene_mat_full[matched_projected_genes, , drop = FALSE]

# ------------------------------------------------------------
# 5) Align test genes to model genes
# ------------------------------------------------------------

test_aligned <- matrix(
  NA_real_,
  nrow = length(genes_used),
  ncol = ncol(test_gene_mat),
  dimnames = list(genes_used, colnames(test_gene_mat))
)

matched_test_genes <- intersect(genes_used, rownames(test_gene_mat))
missing_test_genes <- setdiff(genes_used, rownames(test_gene_mat))

if (length(matched_test_genes) > 0) {
  test_aligned[matched_test_genes, ] <- test_gene_mat[matched_test_genes, , drop = FALSE]
}

# ------------------------------------------------------------
# 6) Impute missing values using training medians
# ------------------------------------------------------------
for (g in genes_used) {
  idx_na <- is.na(test_aligned[g, ])
  if (any(idx_na)) {
    test_aligned[g, idx_na] <- train_gene_medians[g]
  }
}

# Final safety check after imputation
if (anyNA(test_aligned)) {
  stop("NA values remain in aligned test matrix after imputation.")
}

# ------------------------------------------------------------
# 7) Standardize test set using training means and SDs
#    IMPORTANT: use training statistics only
# ------------------------------------------------------------
test_scaled <- sweep(test_aligned, 1, train_gene_means, FUN = "-")
test_scaled <- sweep(test_scaled, 1, train_gene_sds, FUN = "/")
storage.mode(test_scaled) <- "double"

if (any(!is.finite(test_scaled))) {
  bad_n <- sum(!is.finite(test_scaled))
  log_msg("[WARN] Non-finite values detected after scaling: ", bad_n, " -> replacing with 0")
  test_scaled[!is.finite(test_scaled)] <- 0
}

X_test <- t(test_scaled)
storage.mode(X_test) <- "double"

# ------------------------------------------------------------
# 8) Save objects to state
# ------------------------------------------------------------
state$test_meta <- meta_sub
state$test_gene_mat_gxs <- test_aligned              # aligned, imputed, not scaled
state$test_gene_mat_scaled_gxs <- test_scaled        # aligned, imputed, scaled
state$X_test <- X_test                               # sample x gene
state$test_missing_genes <- missing_test_genes
state$test_matched_genes <- matched_test_genes
state$test_region_used <- TEST_REGION_USE

# Optional bookkeeping
state$test_scaling_means <- train_gene_means
state$test_scaling_sds <- train_gene_sds
state$test_impute_medians <- train_gene_medians

# ------------------------------------------------------------
# 9) Save reproducibility tables
# ------------------------------------------------------------
test_tbl_dir <- file.path(CFG$out_dir, "tables")
.dir_create(test_tbl_dir)

# Sample-level metadata + labels
meta_out <- data.table::as.data.table(meta_sub, keep.rownames = "SampleID")
meta_out[, y_test_chr := state$y_test_chr]
meta_out[, y_test := state$y_test]

data.table::fwrite(
  meta_out,
  file.path(test_tbl_dir, paste0(CFG$test_gse_id, "_", TEST_REGION_USE, "_test_metadata.csv"))
)

# Gene matching summary
gene_match_dt <- data.table::data.table(
  GeneSymbol = genes_used,
  matched_in_test = genes_used %in% matched_test_genes,
  impute_median = as.numeric(train_gene_medians[genes_used]),
  train_mean = as.numeric(train_gene_means[genes_used]),
  train_sd = as.numeric(train_gene_sds[genes_used])
)

data.table::fwrite(
  gene_match_dt,
  file.path(test_tbl_dir, paste0(CFG$test_gse_id, "_", TEST_REGION_USE, "_gene_alignment_summary.csv"))
)

# Aligned non-scaled gene x sample matrix
test_aligned_out <- data.table::as.data.table(
  cbind(GeneSymbol = rownames(test_aligned), test_aligned)
)
data.table::fwrite(
  test_aligned_out,
  file.path(test_tbl_dir, paste0(CFG$test_gse_id, "_", TEST_REGION_USE, "_test_aligned_imputed_gene_x_sample.csv"))
)

# Scaled gene x sample matrix 
test_scaled_out <- data.table::as.data.table(
  cbind(GeneSymbol = rownames(test_scaled), test_scaled)
)
data.table::fwrite(
  test_scaled_out,
  file.path(test_tbl_dir, paste0(CFG$test_gse_id, "_", TEST_REGION_USE, "_test_aligned_scaled_gene_x_sample.csv"))
)

# Final model input matrix: sample x gene
X_test_out <- data.table::as.data.table(
  cbind(SampleID = rownames(X_test), X_test)
)
data.table::fwrite(
  X_test_out,
  file.path(test_tbl_dir, paste0(CFG$test_gse_id, "_", TEST_REGION_USE, "_X_test_sample_x_gene.csv"))
)

log_msg("[SUMMARY] External test region: ", TEST_REGION_USE)
log_msg("[SUMMARY] Test samples retained: ", nrow(X_test))
log_msg("[SUMMARY] AD samples: ", sum(state$y_test == 1))
log_msg("[SUMMARY] Control samples: ", sum(state$y_test == 0))
log_msg("[SUMMARY] Model genes: ", length(genes_used))
log_msg("[SUMMARY] Matched test genes: ", length(matched_test_genes))
log_msg("[SUMMARY] Missing test genes imputed from training medians: ", length(missing_test_genes))

log_msg("DONE EXTERNAL TEST PREPARATION: ", CFG$test_gse_id, " [", TEST_REGION_USE, "]")
