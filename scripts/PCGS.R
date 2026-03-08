#!/usr/bin/env Rscript

# PCGS pipeline (GSE118553, Frontal)
# - limma DE for 3 contrasts
# - core = DEG(AsymAD vs Control) ∩ (DEG(AD vs Control) ∪ DEG(AD vs AsymAD))
# - probe -> gene symbol collapse
# - export core expression matrix
# - PCA expansion to expand the feature space

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(matrixStats)
  library(dplyr)
  library(stringr)
  library(tibble)
})

source(file.path("R", "utils.R"))

# ----------------------------
# Parameters (edit here)
# ----------------------------
params <- list(
  gse_id = "GSE118553",
  out_dir = "outputs",
  # phenodata column guesses
  group_col_patterns  = c("disease", "diagnos", "status", "group", "phenotype", "characteristics"),
  region_col_patterns = c("region", "brain", "tissue", "area", "cortex", "pfc", "brodmann"),
  # region matching
  frontal_regex = "frontal|pfc|prefrontal|dlpfc|dorsolateral",
  # DE thresholds
  fdr = 0.05,
  lfc = log2(1.2),
  # expression filter (very light; safe defaults)
  expr_q = 0.25,
  min_prop = 0.20,
  var_q = 0.10,
  # PCA expansion
  pc_cumvar = 0.80,
  cor_grid = seq(0.60, 0.15, by = -0.01)
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Download
# ----------------------------
message("[1/7] Downloading ", params$gse_id)

gset_list <- getGEO(params$gse_id, GSEMatrix = TRUE, getGPL = FALSE)
gset <- gset_list[[1]]

expr <- exprs(gset)   # features x samples
pdat <- pData(gset)

# ----------------------------
# 2) Parse group + region, keep frontal
# ----------------------------
message("[2/7] Parsing phenotype columns")

group_col  <- pick_col(pdat, params$group_col_patterns, required = TRUE)
region_col <- pick_col(pdat, params$region_col_patterns, required = TRUE)

grp <- normalize_group(pdat[[group_col]])
keep <- !is.na(grp)

expr <- expr[, keep, drop = FALSE]
pdat <- pdat[keep, , drop = FALSE]
grp  <- factor(grp[keep], levels = c("Control", "AsymAD", "AD"))

region_text <- tolower(as.character(pdat[[region_col]]))
is_frontal <- str_detect(region_text, params$frontal_regex)

message("Frontal matched samples: ", sum(is_frontal), " / ", length(is_frontal))
message("Group counts (frontal):")
print(table(grp[is_frontal]))

expr_f <- expr[, is_frontal, drop = FALSE]
pdat_f <- pdat[is_frontal, , drop = FALSE]
grp_f  <- droplevels(grp[is_frontal])

if (any(table(grp_f) < 3)) {
  warning("Some groups have <3 samples after frontal filtering; DE may be unstable.")
}

# ----------------------------
# 3) Light feature filtering (optional but recommended)
# ----------------------------
message("[3/7] Filtering low-signal features")
q_ref <- quantile(expr_f, probs = params$expr_q, na.rm = TRUE)
min_n <- ceiling(ncol(expr_f) * params$min_prop)
keep_expr <- rowSums(expr_f > q_ref) >= min_n

v <- rowVars(expr_f, na.rm = TRUE)
keep_var <- v > quantile(v, params$var_q, na.rm = TRUE)

expr_f2 <- expr_f[keep_expr & keep_var, , drop = FALSE]
message("Remaining features: ", nrow(expr_f2))

# ----------------------------
# 4) limma DE
# ----------------------------
message("[4/7] Running limma")

design <- model.matrix(~ 0 + grp_f)
colnames(design) <- make.names(colnames(design))

fit <- lmFit(expr_f2, design)
ct <- makeContrasts(
  AsymAD_vs_Control = grp_fAsymAD - grp_fControl,
  AD_vs_Control     = grp_fAD - grp_fControl,
  AD_vs_AsymAD      = grp_fAD - grp_fAsymAD,
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, ct))

get_deg <- function(fit_obj, coef_name, fdr, lfc) {
  topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("feature") %>%
    as_tibble() %>%
    filter(adj.P.Val < fdr, abs(logFC) >= lfc)
}

deg1 <- get_deg(fit2, "AsymAD_vs_Control", params$fdr, params$lfc)
deg2 <- get_deg(fit2, "AD_vs_Control",     params$fdr, params$lfc)
deg3 <- get_deg(fit2, "AD_vs_AsymAD",      params$fdr, params$lfc)

message("DEG counts: ")
message("  AsymAD vs Control: ", nrow(deg1))
message("  AD vs Control:     ", nrow(deg2))
message("  AD vs AsymAD:      ", nrow(deg3))

core_probes <- intersect(deg1$feature, union(deg2$feature, deg3$feature))
message("Core probes: ", length(core_probes))

# Save DE tables
write.table(deg1, file = file.path(params$out_dir, "deg_AsymAD_vs_Control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg2, file = file.path(params$out_dir, "deg_AD_vs_Control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg3, file = file.path(params$out_dir, "deg_AD_vs_AsymAD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(core_probes, file = file.path(params$out_dir, "core_probes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ----------------------------
# 5) Probe -> gene symbol collapse (GPL)
# ----------------------------
message("[5/7] Collapsing probes to gene symbols")

gpl_id <- annotation(gset)
if (is.na(gpl_id) || gpl_id == "") {
  stop("Could not infer GPL ID from annotation(gset).")
}

gpl <- getGEO(gpl_id, AnnotGPL = TRUE)
gpl_tbl <- Table(gpl)

symbol_col <- pick_symbol_col_from_gpl(gpl_tbl)
if (is.na(symbol_col)) {
  stop("Could not find a gene symbol column in GPL table.")
}

probe_col <- colnames(gpl_tbl)[1]
map_df <- gpl_tbl %>%
  dplyr::select(probe = dplyr::all_of(probe_col), symbol_raw = dplyr::all_of(symbol_col)) %>%
  dplyr::mutate(symbol = clean_symbol(symbol_raw)) %>%
  dplyr::filter(!is.na(symbol))

map_df <- map_df %>% dplyr::filter(probe %in% rownames(expr_f2))
if (nrow(map_df) == 0) {
  stop("No overlap between expression probe IDs and GPL probe IDs.")
}

expr_mapped <- expr_f2[map_df$probe, , drop = FALSE]
expr_gene_all <- limma::avereps(expr_mapped, ID = map_df$symbol) # genes x samples

# Core probes -> core gene symbols
core_genes <- map_df$symbol[match(core_probes, map_df$probe)]
core_genes <- unique(core_genes[!is.na(core_genes)])

write.table(core_genes, file = file.path(params$out_dir, "core_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Export core expression matrix (genes x samples) + group label
expr_core <- expr_gene_all[intersect(core_genes, rownames(expr_gene_all)), , drop = FALSE]

meta_out <- tibble(
  sample = colnames(expr_core),
  group = as.character(grp_f)
)
write.csv(meta_out, file.path(params$out_dir, "sample_metadata_frontal.csv"), row.names = FALSE)

core_mat <- cbind(GeneSymbol = rownames(expr_core), as.data.frame(expr_core))
write.csv(core_mat, file.path(params$out_dir, "core_expression_gene_x_sample.csv"), row.names = FALSE)

# 6) is PCA expansion but it will be in following codes so skipped here.
# ----------------------------
# 7) Save session object
# ----------------------------
message("[7/7] Saving RDS")

saveRDS(list(
  params = params,
  group_col = group_col,
  region_col = region_col,
  frontal_regex = params$frontal_regex,
  expr_frontal_filtered = expr_f2,
  group_frontal = grp_f,
  deg = list(
    AsymAD_vs_Control = deg1,
    AD_vs_Control = deg2,
    AD_vs_AsymAD = deg3
  ),
  core_probes = core_probes,
  core_genes  = core_genes,
  expr_gene_all = expr_gene_all,
  expr_core = expr_core
), file = file.path(params$out_dir, "pcgs_frontal_pipeline.rds"))

message("Done. Outputs written to: ", params$out_dir)
