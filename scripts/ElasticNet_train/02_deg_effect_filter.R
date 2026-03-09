# scripts/02_deg_effect_filter.R
source(file.path("scripts/ElasticNet_train/00_setup.R"))

log_msg("RUNNING FILE: scripts/ElasticNet_train/02_deg_effect_filter.R")

# ------------------------------------------------------------
# 1) Load projected expression matrix
# ------------------------------------------------------------

expr_path <- file.path("data","processed","gse33000_expr_projected_gene_x_sample.csv")

if (!file.exists(expr_path)) {
  stop("Missing projected expression matrix: ", expr_path,
       "\nRun scripts/01_build_expr_gse33000_with_gpl.R first.")
}

log_msg("[1/3] Loading projected expression: ", expr_path)

expr_dt <- data.table::fread(expr_path)

stopifnot("GeneSymbol" %in% names(expr_dt))

expr <- as.data.frame(expr_dt[, -1])
rownames(expr) <- expr_dt$GeneSymbol


# ------------------------------------------------------------
# 2) Load phenotype metadata
# ------------------------------------------------------------

meta_path <- file.path("data","processed","GSE33000_pheno.csv")

if (!file.exists(meta_path)) {
  stop("Missing metadata file: ", meta_path)
}

log_msg("[2/3] Loading metadata: ", meta_path)

meta_raw <- data.table::fread(meta_path)
meta_raw <- as.data.frame(meta_raw)

if (!("Sample" %in% names(meta_raw))) {
  stop("Metadata must contain column 'Sample'")
}

# ------------------------------------------------------------
# 3) Detect AD / Control labels
# ------------------------------------------------------------

log_msg("Detecting AD / Control labels from metadata")

pheno_text <- apply(meta_raw,1,function(x) paste(as.character(x),collapse=" | "))

condition <- ifelse(
  grepl("Alzheimer",pheno_text,ignore.case=TRUE),"AD",
  ifelse(grepl("non[- ]?demented|control",pheno_text,ignore.case=TRUE),"Control",NA)
)

meta <- data.frame(
  sample_id = meta_raw$Sample,
  condition = condition,
  stringsAsFactors = FALSE
)

meta <- meta[!is.na(meta$condition),]

log_msg("Condition counts:")
print(table(meta$condition))


# ------------------------------------------------------------
# 4) Match samples between expression and metadata
# ------------------------------------------------------------

common_samples <- intersect(colnames(expr), meta$sample_id)

if (length(common_samples) < 20) {
  stop("Too few matched samples: ", length(common_samples))
}

expr <- expr[, common_samples, drop=FALSE]

meta <- meta[match(common_samples, meta$sample_id), , drop=FALSE]

stopifnot(identical(colnames(expr), meta$sample_id))

log_msg("Matched samples: ", length(common_samples))
log_msg("AD samples: ", sum(meta$condition=="AD"))
log_msg("Control samples: ", sum(meta$condition=="Control"))


# ------------------------------------------------------------
# 5) Run limma DEG
# ------------------------------------------------------------

log_msg("[3/3] Running limma DEG")

group <- factor(meta$condition, levels=c("Control","AD"))

design <- model.matrix(~group)

fit <- limma::lmFit(expr, design)
fit <- limma::eBayes(fit)

deg <- limma::topTable(
  fit,
  coef=2,
  number=Inf,
  sort.by="P"
)

deg$GeneSymbol <- rownames(deg)

out_deg <- file.path("outputs","tables","deg_gse33000_postproj.csv")

dir.create(dirname(out_deg),recursive=TRUE,showWarnings=FALSE)

data.table::fwrite(deg,out_deg)

log_msg("Saved DEG table: ",out_deg)


# ------------------------------------------------------------
# 6) Effect-size filtering
# ------------------------------------------------------------

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.2

deg_f <- deg[
  deg$adj.P.Val < FDR_CUTOFF &
    abs(deg$logFC) > LOGFC_CUTOFF,
]

out_deg_f <- file.path("outputs","tables","deg_filtered_genes_postproj.csv")

data.table::fwrite(deg_f,out_deg_f)

log_msg("Saved filtered DEG table: ",out_deg_f)


# ------------------------------------------------------------
# 7) Save filtered expression matrix (for LASSO)
# ------------------------------------------------------------

keep <- intersect(rownames(expr),deg_f$GeneSymbol)

expr_f <- expr[keep,,drop=FALSE]

out_expr_f <- file.path(
  "data","processed",
  "gse33000_expr_filtered_gene_x_sample.csv"
)

expr_f_dt <- data.table::as.data.table(expr_f,keep.rownames="GeneSymbol")

data.table::fwrite(expr_f_dt,out_expr_f)

log_msg("Saved filtered expression matrix: ",out_expr_f)

log_msg("Filtered genes: ",nrow(expr_f))
log_msg("DONE")
