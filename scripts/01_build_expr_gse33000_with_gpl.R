log_msg("RUNNING FILE: scripts/01_build_expr_gse33000_with_gpl.R")

if (!exists("state")) state <- new.env(parent = emptyenv())
GSE_ID <- "GSE33000"

OUT_EXPR_PROBE_RDS <- file.path("data", "processed", "GSE33000_expr_probe.rds")
OUT_PHENO_CSV      <- file.path("data", "processed", "GSE33000_pheno.csv")
OUT_MAP_CSV        <- file.path("data", "processed", "GSE33000_probe2gene.csv")
OUT_EXPR_GENE_CSV  <- file.path("data", "processed", "GSE33000_expr_gene.csv")

# ------------------------------------------------------------
# Load precomputed WGCNA-derived core gene set
# ------------------------------------------------------------
core_gene_path <- CFG$core_gene_csv

if (is.null(core_gene_path) || !nzchar(core_gene_path) || !file.exists(core_gene_path)) {
  stop("Core gene file not found: ", core_gene_path)
}

core_dt <- data.table::fread(core_gene_path)

gene_col_candidates <- c(
  "GeneSymbol", "gene_symbol", "SYMBOL", "symbol", "gene", "Gene"
)

gene_col <- gene_col_candidates[gene_col_candidates %in% names(core_dt)][1]
if (is.na(gene_col) || !nzchar(gene_col)) {
  gene_col <- names(core_dt)[1]
}

core_genes <- unique(trimws(as.character(core_dt[[gene_col]])))
core_genes <- core_genes[!is.na(core_genes) & core_genes != ""]
core_genes <- unique(core_genes)

if (length(core_genes) < 20) {
  stop("Too few core genes loaded from: ", core_gene_path)
}

log_msg("[0/4] Loaded precomputed core genes: ", length(core_genes), " from ", core_gene_path)

log_msg("[1/4] Downloading/reading ", GSE_ID, " via getGEO(GSEMatrix=TRUE)")
gse_list <- GEOquery::getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE)
gse <- gse_list[[1]]

expr_mat <- Biobase::exprs(gse)
pheno_dt <- data.table::as.data.table(Biobase::pData(gse), keep.rownames = TRUE)
data.table::setnames(pheno_dt, "rn", "Sample")

pheno_text <- apply(pheno_dt, 1, function(x) paste(as.character(x), collapse = " | "))
hd_idx <- grepl("Huntington's disease", pheno_text, ignore.case = TRUE)
if (any(hd_idx)) {
  hd_samples <- pheno_dt$Sample[hd_idx]
  pheno_dt <- pheno_dt[!Sample %in% hd_samples]
  expr_mat <- expr_mat[, !colnames(expr_mat) %in% hd_samples, drop = FALSE]
}

data.table::fwrite(pheno_dt, OUT_PHENO_CSV)
saveRDS(expr_mat, OUT_EXPR_PROBE_RDS)

log_msg("[2/4] Detected platform: ", annotation(gse), " -> downloading GPL")
gpl_dt <- data.table::as.data.table(
  GEOquery::Table(GEOquery::getGEO(annotation(gse)))
)

# ------------------------------------------------------------
# Helpers for robust GPL column detection and mapping cleanup
# ------------------------------------------------------------

normalize_colname <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

normalize_probe_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[\r\n\t]+", "", x)
  x <- gsub("\\s+", "", x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'+|'+$", "", x)
  x <- sub("\\.0+$", "", x)   # e.g. 10023807248.0 -> 10023807248
  x
}

normalize_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'+|'+$", "", x)
  x
}

pick_col_robust <- function(nms, candidates, required = FALSE, label = "target column") {
  nms <- as.character(nms)
  nms_norm <- normalize_colname(nms)
  cand_norm <- normalize_colname(candidates)
  
  # 1) exact normalized match
  hit_idx <- match(cand_norm, nms_norm, nomatch = 0L)
  hit_idx <- hit_idx[hit_idx > 0L]
  if (length(hit_idx) > 0L) {
    return(nms[hit_idx[1]])
  }
  
  # 2) fallback: substring match
  for (cand in cand_norm) {
    idx <- which(grepl(cand, nms_norm, fixed = TRUE))
    if (length(idx) > 0L) {
      return(nms[idx[1]])
    }
  }
  
  if (required) {
    stop(
      "Could not detect ", label, ".\nAvailable GPL columns:\n",
      paste(sprintf(" - %s", nms), collapse = "\n")
    )
  }
  
  NA_character_
}

clean_symbol_safe <- function(x) {
  x <- normalize_gene_symbol(x)
  
  # Keep only the first symbol when multiple mappings are present
  x <- sub("\\s*(///|//|;|,|\\|).*$", "", x)
  
  # Remove placeholders / unusable entries
  bad <- is.na(x) |
    x == "" |
    x == "---" |
    toupper(x) %in% c("NA", "N/A", "NULL", "UNKNOWN")
  
  x[bad] <- NA_character_
  x
}

# ------------------------------------------------------------
# Detect probe and symbol columns robustly
# ------------------------------------------------------------
log_msg("[2/4] GPL columns: ", paste(names(gpl_dt), collapse = " | "))

probe_candidates <- c(
  "ID", "ID_REF", "ProbeID", "Probe ID", "SPOT_ID", "Reporter ID"
)

symbol_candidates <- c(
  "GENE_SYMBOL", "Gene Symbol", "GeneSymbol", "Symbol", "ORF",
  "GENE ASSIGNMENT", "Gene Assignment"
)

probe_col <- pick_col_robust(
  names(gpl_dt),
  probe_candidates,
  required = FALSE,
  label = "probe column"
)
if (is.na(probe_col)) probe_col <- names(gpl_dt)[1]

symbol_col <- pick_col_robust(
  names(gpl_dt),
  symbol_candidates,
  required = TRUE,
  label = "gene symbol-like column"
)

log_msg("[2/4] Using probe column: ", probe_col)
log_msg("[2/4] Using symbol column: ", symbol_col)

# ------------------------------------------------------------
# Build a clean probe-to-gene mapping table
# Output columns are standardized to:
#   - ProbeID
#   - GeneSymbol
# ------------------------------------------------------------
map_dt <- data.table::data.table(
  ProbeID    = normalize_probe_id(gpl_dt[[probe_col]]),
  GeneSymbol = clean_symbol_safe(gpl_dt[[symbol_col]])
)

# Drop rows with missing / empty mappings
map_dt <- map_dt[
  !is.na(ProbeID) & ProbeID != "" &
    !is.na(GeneSymbol) & GeneSymbol != ""
]

# Remove exact duplicates
map_dt <- unique(map_dt, by = c("ProbeID", "GeneSymbol"))

log_msg("[2/4] Clean mapping rows: ", nrow(map_dt))
log_msg("[2/4] Unique probes mapped: ", data.table::uniqueN(map_dt$ProbeID))
log_msg("[2/4] Unique genes mapped : ", data.table::uniqueN(map_dt$GeneSymbol))

# Optional debug preview
log_msg("[2/4] Mapping preview:")
print(utils::head(map_dt, 10))

# Save the cleaned mapping table
data.table::fwrite(map_dt, OUT_MAP_CSV)

log_msg("[3/4] Mapping & collapsing probes to gene-level expression")
expr_gene_mat <- collapse_to_gene_maxvar(
  expr_mat,
  map_dt,
  keep_genes = NULL
)

storage.mode(expr_gene_mat) <- "double"

expr_gene_dt <- data.table::as.data.table(expr_gene_mat, keep.rownames = TRUE)
data.table::setnames(expr_gene_dt, "rn", "GeneSymbol")

out_expr_gene_full <- file.path("data", "processed", "gse33000_expr_gene_x_sample.csv")
out_expr_proj      <- file.path("data", "processed", "gse33000_expr_projected_gene_x_sample.csv")
out_r2             <- file.path("data", "processed", "GSE33000_core_gene_R2.csv")
out_sel            <- file.path("data", "processed", "GSE33000_selected_core_genes.txt")
out_sel_r2         <- file.path("data", "processed", "GSE33000_selected_core_gene_R2.csv")

data.table::fwrite(expr_gene_dt, out_expr_gene_full)

# ------------------------------------------------------------
# Project training expression matrix onto the precomputed
# WGCNA-derived core-gene feature space
# ------------------------------------------------------------
core_present <- intersect(core_genes, rownames(expr_gene_mat))

if (length(core_present) < 20) {
  stop("Too few core genes present in expr_gene_mat: ", length(core_present))
}

log_msg("[4/4] Core genes present in training matrix: ", length(core_present))

# Samples x core genes
X <- t(expr_gene_mat[core_present, , drop = FALSE])
storage.mode(X) <- "double"

# Fast and robust missing-value handling before PCA:
# remove all-NA / zero-variance genes, then median-impute remaining NAs
non_all_na <- colSums(!is.na(X)) > 0
X <- X[, non_all_na, drop = FALSE]

gene_sds <- apply(X, 2, sd, na.rm = TRUE)
keep_var <- !is.na(gene_sds) & gene_sds > 0
X <- X[, keep_var, drop = FALSE]

if (ncol(X) < 20) {
  stop("Too few usable core genes remained after NA/variance filtering: ", ncol(X))
}

if (anyNA(X)) {
  col_medians <- apply(X, 2, median, na.rm = TRUE)
  for (j in which(colSums(is.na(X)) > 0)) {
    idx <- is.na(X[, j])
    X[idx, j] <- col_medians[j]
  }
}

# PCA on the core-gene space
pca <- prcomp(X, center = TRUE, scale. = TRUE)

ve <- (pca$sdev^2) / sum(pca$sdev^2)
cumve <- cumsum(ve)

# Use enough PCs to explain 80% variance, with a practical cap
K <- which(cumve >= 0.80)[1]
if (is.na(K)) K <- min(ncol(pca$x), 10)
K <- min(K, ncol(pca$x))

scores <- pca$x[, 1:K, drop = FALSE]

# Compute how well each core gene is represented by the retained PC space
r2 <- numeric(ncol(X))
for (j in seq_len(ncol(X))) {
  yj <- X[, j]
  fit <- lm.fit(x = cbind(1, scores), y = yj)
  sse <- sum(fit$residuals^2)
  sst <- sum((yj - mean(yj))^2)
  r2[j] <- if (sst > 0) 1 - sse / sst else NA_real_
}

r2_dt <- data.table::data.table(
  GeneSymbol = colnames(X),
  R2 = r2
)
r2_dt <- r2_dt[!is.na(R2)][order(-R2)]

# Keep a compact projected feature space for downstream ENet
PROJECT_TOP_N <- 600
selected <- r2_dt[1:min(PROJECT_TOP_N, .N), GeneSymbol]

if (length(selected) < 50) {
  stop("Too few projected genes selected: ", length(selected))
}

expr_proj_mat <- expr_gene_mat[selected, , drop = FALSE]
expr_proj_dt <- data.table::as.data.table(expr_proj_mat, keep.rownames = TRUE)
data.table::setnames(expr_proj_dt, "rn", "GeneSymbol")

data.table::fwrite(expr_proj_dt, out_expr_proj)
data.table::fwrite(r2_dt, out_r2)
data.table::fwrite(r2_dt[GeneSymbol %in% selected], out_sel_r2)
writeLines(selected, out_sel)

log_msg(
  "DONE. genes(full)=", nrow(expr_gene_mat),
  " genes(projected)=", nrow(expr_proj_mat),
  " samples=", ncol(expr_gene_mat),
  " PCs_used=", K
)
