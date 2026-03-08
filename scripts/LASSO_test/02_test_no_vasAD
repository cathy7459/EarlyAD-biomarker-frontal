# scripts/00_split_gse122063_by_region_and_diagnosis.R
# source("config/LASSO_test.R")
# source("R/utils/geo.R")

log_msg("=== SPLIT GSE122063 BY REGION / DIAGNOSIS BEFORE PREPROCESSING ===")

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
})

if (!exists("state")) state <- new.env(parent = emptyenv())

if (!exists("TEST_REGION_COL")) {
  TEST_REGION_COL <- "brain region:ch1"
}

out_dir_raw <- file.path(CFG$out_dir, "split_raw")
dir.create(out_dir_raw, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1) Load GEO
# ------------------------------------------------------------
log_msg("[1/6] Downloading/loading GEO: ", CFG$test_gse_id)

gse_obj <- GEOquery::getGEO(CFG$test_gse_id, GSEMatrix = TRUE)

if (length(gse_obj) > 1) {
  log_msg("[Split] multiple ExpressionSets found: ", length(gse_obj), " -> using first")
  eset <- gse_obj[[1]]
} else {
  eset <- gse_obj[[1]]
}

expr_probe <- Biobase::exprs(eset)   # probes x samples
ph <- Biobase::pData(eset)

log_msg("[Split] raw probe matrix: ", nrow(expr_probe), " probes x ", ncol(expr_probe), " samples")

# ------------------------------------------------------------
# 2) Check phenotype columns
# ------------------------------------------------------------
if (!(CFG$test_diag_col %in% colnames(ph))) {
  hit <- colnames(ph)[tolower(colnames(ph)) == tolower(CFG$test_diag_col)]
  if (length(hit) == 1) {
    CFG$test_diag_col <- hit
  } else {
    stop(
      "Diagnosis column not found. Wanted: ", CFG$test_diag_col,
      " | Available: ", paste(colnames(ph), collapse = ", ")
    )
  }
}

if (!(TEST_REGION_COL %in% colnames(ph))) {
  hit_region <- colnames(ph)[tolower(colnames(ph)) == tolower(TEST_REGION_COL)]
  if (length(hit_region) == 1) {
    TEST_REGION_COL <- hit_region
  } else {
    stop(
      "Region column not found. Wanted: ", TEST_REGION_COL,
      " | Available: ", paste(colnames(ph), collapse = ", ")
    )
  }
}

diag_raw <- as.character(ph[[CFG$test_diag_col]])
diag_lc <- trimws(tolower(diag_raw))

region_raw <- as.character(ph[[TEST_REGION_COL]])
region_lc <- trimws(tolower(region_raw))

# ------------------------------------------------------------
# 3) Define labels BEFORE any preprocessing
# ------------------------------------------------------------
# exact labels in GSE122063:
#   - "Alzheimer's disease"
#   - "Control"
#   - "Vascular dementia"

is_vascular <- diag_lc == "vascular dementia"
is_ad <- diag_lc == "alzheimer's disease"
is_ctrl <- diag_lc == "control"

is_frontal <- grepl("frontal", region_lc)
is_temporal <- grepl("temporal", region_lc)

# split A: frontal + non-vascular AD/Control
keep_frontal_nonvascular <- which(is_frontal & (is_ad | is_ctrl) & !is_vascular)

# split B: temporal + non-vascular AD/Control
keep_temporal_nonvascular <- which(is_temporal & (is_ad | is_ctrl) & !is_vascular)

log_msg("[3/6] Phenotype summary:")
print(table(region_lc, diag_lc, useNA = "ifany"))

log_msg("[3/6] frontal (exclude vascular, keep AD/Control): ", length(keep_frontal_nonvascular))
log_msg("[3/6] temporal (exclude vascular, keep AD/Control): ", length(keep_temporal_nonvascular))

if (length(keep_frontal_nonvascular) == 0) {
  stop("No samples found for frontal + AD/Control + exclude vascular dementia.")
}
if (length(keep_temporal_nonvascular) == 0) {
  stop("No samples found for temporal + AD/Control + exclude vascular dementia.")
}

# ------------------------------------------------------------
# 4) Subset probe-level matrices FIRST
# ------------------------------------------------------------
expr_probe_frontal <- expr_probe[, keep_frontal_nonvascular, drop = FALSE]
expr_probe_temporal <- expr_probe[, keep_temporal_nonvascular, drop = FALSE]

meta_frontal <- ph[keep_frontal_nonvascular, , drop = FALSE]
meta_temporal <- ph[keep_temporal_nonvascular, , drop = FALSE]

# save metadata in state too
state$gse122063_frontal_meta <- meta_frontal
state$gse122063_temporal_meta <- meta_temporal

# save raw probe-level subsets
fwrite(
  data.table(
    ProbeID = rownames(expr_probe_frontal),
    as.data.table(expr_probe_frontal)
  ),
  file.path(out_dir_raw, paste0(CFG$test_gse_id, "_frontal_AD_Control_excluding_vascular_probe_x_sample.csv"))
)

fwrite(
  data.table(
    ProbeID = rownames(expr_probe_temporal),
    as.data.table(expr_probe_temporal)
  ),
  file.path(out_dir_raw, paste0(CFG$test_gse_id, "_temporal_AD_Control_excluding_vascular_probe_x_sample.csv"))
)

fwrite(
  as.data.table(meta_frontal, keep.rownames = "sample_id"),
  file.path(out_dir_raw, paste0(CFG$test_gse_id, "_frontal_AD_Control_excluding_vascular_metadata.csv"))
)

fwrite(
  as.data.table(meta_temporal, keep.rownames = "sample_id"),
  file.path(out_dir_raw, paste0(CFG$test_gse_id, "_temporal_AD_Control_excluding_vascular_metadata.csv"))
)

log_msg("[4/6] Saved probe-level split matrices and metadata")

# ------------------------------------------------------------
# 5) Probe -> gene collapse AFTER split
# ------------------------------------------------------------
platform_id <- annotation(eset)
if (is.null(platform_id) || platform_id == "") {
  platform_id <- tryCatch(as.character(Biobase::experimentData(eset)@annotation),
                          error = function(e) "")
}
if (is.null(platform_id) || platform_id == "") {
  stop("Cannot determine platform ID.")
}

log_msg("[5/6] platform detected: ", platform_id)

test_anno <- get_probe_gene_mapping(platform_id)

collapsed_frontal <- collapse_to_gene_maxvar(
  expr_probe_frontal,
  test_anno,
  keep_genes = NULL
)

collapsed_temporal <- collapse_to_gene_maxvar(
  expr_probe_temporal,
  test_anno,
  keep_genes = NULL
)

expr_gene_frontal_gxs <- collapsed_frontal$expr_gene
expr_gene_temporal_gxs <- collapsed_temporal$expr_gene

storage.mode(expr_gene_frontal_gxs) <- "double"
storage.mode(expr_gene_temporal_gxs) <- "double"

log_msg("[5/6] frontal collapsed gene matrix: ",
        nrow(expr_gene_frontal_gxs), " genes x ", ncol(expr_gene_frontal_gxs), " samples")

log_msg("[5/6] temporal collapsed gene matrix: ",
        nrow(expr_gene_temporal_gxs), " genes x ", ncol(expr_gene_temporal_gxs), " samples")

# ------------------------------------------------------------
# 6) Save sample x gene matrices
# ------------------------------------------------------------
expr_gene_frontal_sxg <- t(expr_gene_frontal_gxs)
expr_gene_temporal_sxg <- t(expr_gene_temporal_gxs)

storage.mode(expr_gene_frontal_sxg) <- "double"
storage.mode(expr_gene_temporal_sxg) <- "double"

frontal_dt <- data.table(
  sample_id = rownames(expr_gene_frontal_sxg),
  as.data.table(expr_gene_frontal_sxg)
)

temporal_dt <- data.table(
  sample_id = rownames(expr_gene_temporal_sxg),
  as.data.table(expr_gene_temporal_sxg)
)

out_frontal <- file.path(
  out_dir_raw,
  paste0(CFG$test_gse_id, "_frontal_AD_Control_excluding_vascular_sample_x_gene.csv")
)

out_temporal <- file.path(
  out_dir_raw,
  paste0(CFG$test_gse_id, "_temporal_AD_Control_excluding_vascular_sample_x_gene.csv")
)

fwrite(frontal_dt, out_frontal)
fwrite(temporal_dt, out_temporal)

log_msg("[6/6] Saved sample x gene matrix: ", out_frontal)
log_msg("[6/6] Saved sample x gene matrix: ", out_temporal)

# ------------------------------------------------------------
# 7) Store in state for downstream use
# ------------------------------------------------------------
state$gse122063_frontal_samples <- colnames(expr_probe_frontal)
state$gse122063_temporal_samples <- colnames(expr_probe_temporal)

state$gse122063_frontal_expr_probe <- expr_probe_frontal
state$gse122063_temporal_expr_probe <- expr_probe_temporal

state$gse122063_frontal_expr_gene_gxs <- expr_gene_frontal_gxs
state$gse122063_temporal_expr_gene_gxs <- expr_gene_temporal_gxs

state$gse122063_frontal_expr_gene_sxg <- expr_gene_frontal_sxg
state$gse122063_temporal_expr_gene_sxg <- expr_gene_temporal_sxg

log_msg("DONE: split before preprocessing completed.")
