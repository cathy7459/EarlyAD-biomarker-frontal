log_msg("RUNNING FILE: scripts/02_deg_effect_filter.R")

expr_path <- file.path("data", "processed", "gse33000_expr_projected_gene_x_sample.csv")
if (!file.exists(expr_path)) {
  stop("Missing projected expression matrix: ", expr_path,
       "\nRun scripts/01_build_expr_gse33000_with_gpl.R first.")
}

expr_dt <- data.table::fread(expr_path)
stopifnot("GeneSymbol" %in% names(expr_dt))
expr <- as.data.frame(expr_dt[, -1])
rownames(expr) <- expr_dt$GeneSymbol
expr[] <- lapply(expr, as.numeric)

meta_path <- file.path("data", "processed", "GSE33000_pheno.csv")
if (!file.exists(meta_path)) stop("Missing metadata file: ", meta_path)
meta_raw <- data.table::fread(meta_path)
meta_raw <- as.data.frame(meta_raw)
if (!("Sample" %in% names(meta_raw))) stop("Metadata must contain column 'Sample'")

pheno_text <- apply(meta_raw, 1, function(x) paste(as.character(x), collapse = " | "))
condition <- ifelse(
  grepl("Alzheimer", pheno_text, ignore.case = TRUE), "AD",
  ifelse(grepl("non[- ]?demented|control", pheno_text, ignore.case = TRUE), "Control", NA)
)
meta <- data.frame(sample_id = meta_raw$Sample, condition = condition, stringsAsFactors = FALSE)
meta <- meta[!is.na(meta$condition), ]

common_samples <- intersect(colnames(expr), meta$sample_id)
if (length(common_samples) < 20) stop("Too few matched samples: ", length(common_samples))
expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
stopifnot(identical(colnames(expr), meta$sample_id))

group <- factor(meta$condition, levels = c("Control", "AD"))
design <- model.matrix(~group)
fit <- limma::lmFit(expr, design)
fit <- limma::eBayes(fit)
deg <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg$GeneSymbol <- rownames(deg)

out_deg <- file.path("outputs", "tables", "deg_gse33000_postproj.csv")
out_deg_f <- file.path("outputs", "tables", "deg_filtered_genes_postproj.csv")
out_expr_f <- file.path("data", "processed", "gse33000_expr_filtered_gene_x_sample.csv")
.dir_create(dirname(out_deg))

data.table::fwrite(deg, out_deg)

deg_f <- deg[deg$adj.P.Val < FDR_CUTOFF & abs(deg$logFC) > LOGFC_CUTOFF, ]
data.table::fwrite(deg_f, out_deg_f)

keep <- intersect(rownames(expr), deg_f$GeneSymbol)
expr_f <- expr[keep, , drop = FALSE]
expr_f_dt <- data.table::as.data.table(expr_f, keep.rownames = "GeneSymbol")
data.table::fwrite(expr_f_dt, out_expr_f)

log_msg("Saved DEG table: ", out_deg)
log_msg("Saved filtered DEG table: ", out_deg_f)
log_msg("Saved filtered expression matrix: ", out_expr_f)
log_msg("Filtered genes: ", nrow(expr_f))
log_msg("DONE")
