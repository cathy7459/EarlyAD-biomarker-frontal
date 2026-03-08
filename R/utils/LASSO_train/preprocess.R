suppressPackageStartupMessages({
  library(data.table)
})

collapse_probes_by_maxvar <- function(dt, gene_col, sample_cols) {
  # dt must include gene_col and sample_cols
  dt[, probe_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]
  data.table::setorder(dt, get(gene_col), -probe_var)
  dt_best <- dt[, .SD[1], by = gene_col]
  dt_best[]
}

to_gene_x_sample <- function(dt_best, gene_col, sample_cols) {
  mat <- as.data.frame(dt_best[, ..sample_cols])
  rownames(mat) <- dt_best[[gene_col]]
  mat
}
