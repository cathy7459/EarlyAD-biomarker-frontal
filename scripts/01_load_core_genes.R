if (!exists("state")) state <- new.env(parent = emptyenv())
state$core_genes <- NULL

core_gene_path <- NULL

if (exists("CFG") && !is.null(CFG$core_gene_csv) && nzchar(CFG$core_gene_csv)) {
  if (file.exists(CFG$core_gene_csv)) {
    core_gene_path <- CFG$core_gene_csv
  }
}

if (is.null(core_gene_path)) {
  candidate_paths <- c(
    file.path("data", "processed", "WGCNA_significant_genes.csv"),
    file.path("outputs", "data", "processed", "WGCNA_significant_genes.csv"),
    file.path("outputs", "tables", "WGCNA_significant_genes.csv")
  )
  hit <- candidate_paths[file.exists(candidate_paths)][1]
  if (!is.na(hit)) core_gene_path <- hit
}

if (!is.null(core_gene_path) && file.exists(core_gene_path)) {
  core_dt <- data.table::fread(core_gene_path)
  gene_col <- names(core_dt)[match(
    TRUE,
    tolower(names(core_dt)) %in% c("genesymbol", "gene", "symbol", "hgnc_symbol", "gene_symbol")
  )]
  if (is.na(gene_col)) gene_col <- names(core_dt)[1]

  core_genes <- unique(trimws(as.character(core_dt[[gene_col]])))
  core_genes <- core_genes[!is.na(core_genes) & core_genes != ""]

  state$core_genes <- core_genes
  CFG$core_gene_csv <- core_gene_path

  log_msg("[Core] loaded ", length(core_genes),
          " genes from ", core_gene_path,
          " (col=", gene_col, ")")
} else {
  log_msg("[Core] core gene file not found -> proceeding WITHOUT core filtering (not recommended).")
}
