suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

run_limma_deg <- function(expr_gene_x_sample, meta, condition_col = "condition") {
  # expr: genes x samples numeric matrix/data.frame
  # meta: data.frame with sample_id matching colnames(expr) and condition (AD/Control)

  stopifnot(all(c("sample_id", condition_col) %in% colnames(meta)))

  # align
  meta <- meta[match(colnames(expr_gene_x_sample), meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id)) stop("Sample metadata does not cover all expression columns.")
  group <- factor(meta[[condition_col]])
  if (!all(levels(group) %in% c("AD", "Control"))) {
    # allow any ordering but require both
    if (length(unique(group)) != 2) stop("Need exactly 2 conditions for limma.")
  }

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(expr_gene_x_sample, design)
  # Contrast: AD - Control (if names differ, take 2nd - 1st in alphabetical order)
  lev <- colnames(design)
  if (all(c("AD", "Control") %in% lev)) {
    contr <- makeContrasts(ADvsControl = AD - Control, levels = design)
  } else {
    contr <- makeContrasts(contrast = lev[2] - lev[1], levels = design)
    colnames(contr) <- "Group2vsGroup1"
  }
  fit2 <- eBayes(contrasts.fit(fit, contr))

  tt <- topTable(fit2, number = Inf, sort.by = "P")
  dt <- data.table::as.data.table(tt, keep.rownames = "GeneSymbol")
  dt
}
