log_msg("=== SPLIT GSE122063 BY REGION / DIAGNOSIS BEFORE PREPROCESSING ===")

if (!exists("state")) state <- new.env(parent = emptyenv())
TEST_REGION_COL <- "brain region:ch1"
out_dir_raw <- file.path(CFG$out_dir, "split_raw")
.dir_create(out_dir_raw)

log_msg("[1/2] Downloading/loading GEO: ", CFG$test_gse_id)
gse_obj <- GEOquery::getGEO(CFG$test_gse_id, GSEMatrix = TRUE)
eset <- gse_obj[[1]]
expr_probe <- Biobase::exprs(eset)
ph <- Biobase::pData(eset)

if (!(CFG$test_diag_col %in% colnames(ph))) {
  hit <- colnames(ph)[tolower(colnames(ph)) == tolower(CFG$test_diag_col)]
  if (length(hit) == 1) CFG$test_diag_col <- hit else stop("Diagnosis column not found.")
}
if (!(TEST_REGION_COL %in% colnames(ph))) {
  hit_region <- colnames(ph)[tolower(colnames(ph)) == tolower(TEST_REGION_COL)]
  if (length(hit_region) == 1) TEST_REGION_COL <- hit_region else stop("Region column not found.")
}

diag_lc <- trimws(tolower(as.character(ph[[CFG$test_diag_col]])))
region_lc <- trimws(tolower(as.character(ph[[TEST_REGION_COL]])))

is_vascular <- diag_lc == "vascular dementia"
is_ad <- diag_lc == "alzheimer's disease"
is_ctrl <- diag_lc == "control"
is_frontal <- grepl("frontal", region_lc)
is_temporal <- grepl("temporal", region_lc)

keep_frontal_nonvascular <- which(is_frontal & (is_ad | is_ctrl) & !is_vascular)
keep_temporal_nonvascular <- which(is_temporal & (is_ad | is_ctrl) & !is_vascular)

expr_probe_frontal <- expr_probe[, keep_frontal_nonvascular, drop = FALSE]
expr_probe_temporal <- expr_probe[, keep_temporal_nonvascular, drop = FALSE]
meta_frontal <- ph[keep_frontal_nonvascular, , drop = FALSE]
meta_temporal <- ph[keep_temporal_nonvascular, , drop = FALSE]

platform_id <- annotation(eset)
test_anno_raw <- get_probe_gene_mapping(platform_id)

# ------------------------------------------------------------
# Normalize mapping table to the format required by
# collapse_to_gene_maxvar(): ProbeID, GeneSymbol
# ------------------------------------------------------------
normalize_probe_gene_map <- function(map_obj) {
  map_dt <- data.table::as.data.table(map_obj)
  
  if (nrow(map_dt) == 0) {
    stop("Probe-gene mapping table is empty.")
  }
  
  original_names <- names(map_dt)
  normalized_names <- tolower(trimws(original_names))
  
  # Candidate column names for probe IDs
  probe_candidates <- c(
    "probeid", "probe_id", "id", "id_ref", "probe", "probesetid",
    "reporter id", "reporter_id", "name"
  )
  
  # Candidate column names for gene symbols
  gene_candidates <- c(
    "genesymbol", "gene_symbol", "symbol", "gene symbol",
    "gene assignment", "gene", "symbols"
  )
  
  probe_hit <- original_names[match(probe_candidates, normalized_names, nomatch = 0)]
  probe_hit <- probe_hit[nzchar(probe_hit)][1]
  
  gene_hit <- original_names[match(gene_candidates, normalized_names, nomatch = 0)]
  gene_hit <- gene_hit[nzchar(gene_hit)][1]
  
  if (is.na(probe_hit) || !nzchar(probe_hit)) {
    stop(
      "Could not identify probe ID column in mapping table.\n",
      "Available columns: ", paste(original_names, collapse = ", ")
    )
  }
  
  if (is.na(gene_hit) || !nzchar(gene_hit)) {
    stop(
      "Could not identify gene symbol column in mapping table.\n",
      "Available columns: ", paste(original_names, collapse = ", ")
    )
  }
  
  map_dt <- map_dt[, .(
    ProbeID = as.character(get(probe_hit)),
    GeneSymbol = as.character(get(gene_hit))
  )]
  
  # Clean strings
  map_dt[, ProbeID := trimws(ProbeID)]
  map_dt[, GeneSymbol := trimws(GeneSymbol)]
  
  # Remove obvious empties
  map_dt <- map_dt[
    !is.na(ProbeID) & nzchar(ProbeID) &
      !is.na(GeneSymbol) & nzchar(GeneSymbol)
  ]
  
  # If gene symbol field contains delimiters, keep the first symbol
  map_dt[, GeneSymbol := sub(" ///.*$", "", GeneSymbol)]
  map_dt[, GeneSymbol := sub(" //.*$", "", GeneSymbol)]
  map_dt[, GeneSymbol := sub(";.*$", "", GeneSymbol)]
  map_dt[, GeneSymbol := sub(",.*$", "", GeneSymbol)]
  map_dt[, GeneSymbol := trimws(GeneSymbol)]
  
  map_dt <- map_dt[
    !is.na(GeneSymbol) &
      nzchar(GeneSymbol) &
      GeneSymbol != "---"
  ]
  
  # Keep unique probe-gene pairs
  map_dt <- unique(map_dt, by = c("ProbeID", "GeneSymbol"))
  
  if (!all(c("ProbeID", "GeneSymbol") %in% names(map_dt))) {
    stop("Failed to normalize mapping table to required columns: ProbeID, GeneSymbol")
  }
  
  if (nrow(map_dt) == 0) {
    stop("Normalized probe-gene mapping table is empty after cleaning.")
  }
  
  map_dt
}

test_anno <- normalize_probe_gene_map(test_anno_raw)

log_msg("[2/2] Mapping table normalized: ", nrow(test_anno), " probe-gene pairs")
log_msg("[2/2] Mapping columns: ", paste(names(test_anno), collapse = ", "))

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

# ------------------------------------------------------------
# Robust extractor for gene-level expression output
# Supports either:
#   1) a list containing $expr_gene
#   2) a matrix/data.frame returned directly
# ------------------------------------------------------------
extract_expr_gene_matrix <- function(x, obj_name = "collapsed_object") {
  if (is.list(x) && !is.null(x$expr_gene)) {
    expr_gene <- x$expr_gene
  } else if (is.matrix(x) || is.data.frame(x)) {
    expr_gene <- x
  } else {
    stop(
      obj_name, " is not in a supported format.\n",
      "Expected either a list with $expr_gene or a matrix/data.frame directly."
    )
  }
  
  expr_gene <- as.matrix(expr_gene)
  storage.mode(expr_gene) <- "double"
  
  if (is.null(rownames(expr_gene)) || nrow(expr_gene) == 0) {
    stop(obj_name, ": extracted gene matrix has no valid rownames.")
  }
  if (is.null(colnames(expr_gene)) || ncol(expr_gene) == 0) {
    stop(obj_name, ": extracted gene matrix has no valid colnames.")
  }
  
  expr_gene
}

state$gse122063_frontal_meta <- meta_frontal
state$gse122063_temporal_meta <- meta_temporal
state$gse122063_frontal_expr_probe <- expr_probe_frontal
state$gse122063_temporal_expr_probe <- expr_probe_temporal
state$gse122063_frontal_expr_gene_gxs <- extract_expr_gene_matrix(
  collapsed_frontal,
  "collapsed_frontal"
)
state$gse122063_temporal_expr_gene_gxs <- extract_expr_gene_matrix(
  collapsed_temporal,
  "collapsed_temporal"
)

split_dir <- file.path(CFG$out_dir, "split_processed")
.dir_create(split_dir)

# ------------------------------------------------------------
# Frontal region
# ------------------------------------------------------------

data.table::fwrite(
  as.data.table(state$gse122063_frontal_meta),
  file.path(split_dir, "GSE122063_frontal_metadata.csv")
)

data.table::fwrite(
  data.table::as.data.table(
    cbind(GeneSymbol = rownames(state$gse122063_frontal_expr_gene_gxs),
          state$gse122063_frontal_expr_gene_gxs)
  ),
  file.path(split_dir, "GSE122063_frontal_expr_gene_gxs.csv")
)

data.table::fwrite(
  data.table::as.data.table(
    cbind(ProbeID = rownames(state$gse122063_frontal_expr_probe),
          state$gse122063_frontal_expr_probe)
  ),
  file.path(split_dir, "GSE122063_frontal_expr_probe_x_sample.csv")
)

# ------------------------------------------------------------
# Temporal region
# ------------------------------------------------------------

data.table::fwrite(
  as.data.table(state$gse122063_temporal_meta),
  file.path(split_dir, "GSE122063_temporal_metadata.csv")
)

data.table::fwrite(
  data.table::as.data.table(
    cbind(GeneSymbol = rownames(state$gse122063_temporal_expr_gene_gxs),
          state$gse122063_temporal_expr_gene_gxs)
  ),
  file.path(split_dir, "GSE122063_temporal_expr_gene_gxs.csv")
)

data.table::fwrite(
  data.table::as.data.table(
    cbind(ProbeID = rownames(state$gse122063_temporal_expr_probe),
          state$gse122063_temporal_expr_probe)
  ),
  file.path(split_dir, "GSE122063_temporal_expr_probe_x_sample.csv")
)

log_msg("Saved: frontal + temporal split datasets")

log_msg("DONE: split before preprocessing completed.")
