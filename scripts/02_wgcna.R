#wgcna
#!/usr/bin/env Rscript
# 02_wgcna.R
# ------------------------------------------------------------------------------
# Purpose
#   Run WGCNA on an (expanded) gene x sample expression matrix and export:
#     - Gene-to-module mapping
#     - Module eigengenes
#     - Significant (hub) genes by kME (and optionally trait-associated hubs)
#     - A clean triangle dendrogram with a module color bar
#     - A bundled RDS object with key results
#
# Expected inputs (must exist BEFORE running this script)
#   - expanded_expr : gene x sample numeric matrix (rownames=gene symbols, colnames=sample IDs)
#   - group_f       : (optional) factor/character vector of sample labels (names=sample IDs recommended)
#
# Outputs (written under outputs/)
#   - WGCNA_GeneModuleMap.csv
#   - WGCNA_ModuleEigengenes.csv
#   - WGCNA_significant_genes.csv
#   - WGCNA_ME_trait_correlation.tsv (optional)
#   - WGCNA_ME_trait_correlation_pvalues.tsv (optional)
#   - WGCNA_triangle_with_colorbar.pdf
#   - WGCNA_results.rds
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(WGCNA)
  library(dendextend)
  library(readr)
})
conflicts(detail = TRUE)

# ------------------------------- Helpers --------------------------------------

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

require_objects <- function(obj_names, env = parent.frame()) {
  missing <- obj_names[!vapply(obj_names, exists, logical(1), envir = env)]
  if (length(missing) > 0) {
    stop("Missing required object(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

as_numeric_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("Input must be a matrix or data.frame convertible to a matrix.", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

# ------------------------------- Setup ----------------------------------------

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------- Inputs ---------------------------------------

log_step("WGCNA | input check")
require_objects(c("expanded_expr"))

# Use expanded expression matrix as expr_gene (gene x sample)
expr_gene <- expanded_expr

# Normalize input shape if a GeneSymbol column is present
if (is.data.frame(expr_gene) && "GeneSymbol" %in% colnames(expr_gene)) {
  rownames(expr_gene) <- expr_gene$GeneSymbol
  expr_gene$GeneSymbol <- NULL
  expr_gene <- as.matrix(expr_gene)
}

expr_gene <- as_numeric_matrix(expr_gene)

if (is.null(rownames(expr_gene)) || anyNA(rownames(expr_gene))) {
  stop("expr_gene must have valid rownames (gene symbols).", call. = FALSE)
}
if (is.null(colnames(expr_gene)) || anyNA(colnames(expr_gene))) {
  stop("expr_gene must have valid colnames (sample IDs).", call. = FALSE)
}

# WGCNA expects samples in rows and genes in columns
datExpr <- t(expr_gene)
cat("datExpr dim (samples x genes): ", paste(dim(datExpr), collapse = " x "), "\n")

# Remove problematic samples/genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing genes (bad): ", paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", "), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing samples (bad): ", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", "), "\n")
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# --------------------- 1) Soft-threshold power selection ----------------------
# (robust, avoids power=1)
log_step("WGCNA | pick soft-threshold power")

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

fitTbl <- sft$fitIndices
# WGCNA fitIndices columns typically include:
#   Power, SFT.R.sq, slope, truncated.R.sq, mean.k., median.k., max.k.
# The exact names can vary slightly by version, but SFT.R.sq and mean.k. are common.

# ---- Tuning knobs (edit these if needed)
target_R2 <- 0.85          # stricter than 0.70; common paper choice is 0.80~0.90
min_power <- 6             # signed networks often use >= 6 (prevents power=1 selection)
min_mean_k <- 5            # prevent selecting a power that destroys connectivity

# Candidate powers meeting scale-free fit AND not too low connectivity
cand <- fitTbl[
  fitTbl$SFT.R.sq >= target_R2 &
    fitTbl$Power >= min_power &
    fitTbl$mean.k. >= min_mean_k,
  , drop = FALSE
]

if (nrow(cand) > 0) {
  # Choose the smallest power that satisfies all criteria
  softPower <- min(cand$Power)
} else {
  # Fallback strategy:
  #  1) choose the best R2 among powers >= min_power
  sub <- fitTbl[fitTbl$Power >= min_power, , drop = FALSE]
  if (nrow(sub) == 0) sub <- fitTbl
  
  best <- sub[which.max(sub$SFT.R.sq), , drop = FALSE]
  softPower <- best$Power
  
  # If connectivity is extremely low, back off slightly (optional safety)
  if (!is.na(best$mean.k.) && best$mean.k. < min_mean_k) {
    # pick the power >= min_power with the highest mean connectivity while keeping decent R2
    sub2 <- sub[sub$SFT.R.sq >= 0.70, , drop = FALSE]
    if (nrow(sub2) > 0) {
      softPower <- sub2$Power[which.max(sub2$mean.k.)]
    }
  }
}

cat("Chosen softPower =", softPower, "\n")
# -------------------------- 2) Module detection -------------------------------

log_step("WGCNA | blockwiseModules")

# Parameters are tuned to reduce the number of genes assigned to the "grey" module:
#   - smaller minModuleSize
#   - higher deepSplit
#   - allow PAM reassignments
#   - conservative merging
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  
  minModuleSize = 50,
  deepSplit = 2,
  
  pamRespectsDendro = FALSE,
  reassignThreshold = 0,
  
  mergeCutHeight = 0.30,
  
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

cat("Module counts:\n")
print(table(moduleColors))

log_step("WGCNA | network construction done")

# -------------------------- 3) Module eigengenes ------------------------------

log_step("WGCNA | module eigengenes")
MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)
cat("MEs dim:", nrow(MEs), "samples x", ncol(MEs), "modules\n")

# ----------- 4) Significant genes: hub genes by kME (optional trait) ----------

log_step("WGCNA | extract significant genes (module hubs)")

# kME: correlation between each gene and each module eigengene
kME_mat <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
colnames(kME_mat) <- paste0("kME_", gsub("^ME", "", colnames(MEs)))
kME_mat$GeneSymbol <- rownames(kME_mat)

# Per-gene module assignment
gene_info <- data.frame(
  GeneSymbol  = colnames(datExpr),
  ModuleColor = moduleColors,
  ModuleLabel = moduleLabels,
  stringsAsFactors = FALSE
)

# For each gene, extract kME for its own module
get_own_kME <- function(gene, mod_color) {
  coln <- paste0("kME_", mod_color)
  if (!coln %in% colnames(kME_mat)) return(NA_real_)
  kME_mat[kME_mat$GeneSymbol == gene, coln, drop = TRUE]
}

gene_info$kME_own <- mapply(get_own_kME, gene_info$GeneSymbol, gene_info$ModuleColor)

# Hub gene threshold (adjust as needed)
kME_cut <- 0.85
hub_genes <- gene_info[!is.na(gene_info$kME_own) & abs(gene_info$kME_own) >= kME_cut, ]

sig_genes  <- hub_genes
sig_reason <- "kME_hubs_only"

# Optional: trait-associated module filtering if group_f is provided
p_cut <- 0.01 # BH-FDR cutoff for module-trait association

if (exists("group_f")) {
  log_step("WGCNA | trait-aware selection enabled (group_f found)")
  
  # Align group labels to samples if group_f has names
  if (!is.null(names(group_f))) {
    common <- intersect(rownames(datExpr), names(group_f))
    if (length(common) >= 5) {
      MEs2   <- MEs[common, , drop = FALSE]
      group2 <- group_f[common]
    } else {
      warning("group_f names do not match datExpr rownames well; skipping trait filtering.")
      MEs2 <- NULL
    }
  } else {
    MEs2   <- MEs
    group2 <- group_f
  }
  
  if (!is.null(MEs2)) {
    # Build a simple one-vs-rest trait design (supports Control/AsymAD/AD if present)
    trait <- data.frame(
      Control = as.numeric(group2 == "Control"),
      AsymAD  = as.numeric(group2 == "AsymAD"),
      AD      = as.numeric(group2 == "AD")
    )
    rownames(trait) <- rownames(MEs2)
    
    modTraitCor <- cor(MEs2, trait, use = "pairwise.complete.obs")
    modTraitP   <- corPvalueStudent(modTraitCor, nSamples = nrow(MEs2))
    modTraitFDR <- apply(modTraitP, 2, p.adjust, method = "BH")
    
    # Select modules where ANY trait has BH-FDR < p_cut
    sig_modules <- rownames(modTraitFDR)[apply(modTraitFDR, 1, function(q) any(q < p_cut))]
    sig_modules <- gsub("^ME", "", sig_modules)
    
    if (length(sig_modules) > 0) {
      sig_genes  <- hub_genes[hub_genes$ModuleColor %in% sig_modules, ]
      sig_reason <- paste0("trait_modules_FDR<", p_cut, "_AND_abs(kME)>=", kME_cut)
    } else {
      sig_reason <- paste0("no_trait_modules_FDR<", p_cut, "; used_kME_hubs_only")
    }
  }
}

sig_genes$SelectionReason <- sig_reason

write_csv(sig_genes, file.path(out_dir, "WGCNA_significant_genes.csv"))
cat("Saved significant genes:", nrow(sig_genes), "to outputs/WGCNA_significant_genes.csv\n")

# --------------------- 5) Optional: export trait correlations -----------------

if (exists("group_f")) {
  log_step("WGCNA | export module-trait correlation (optional)")
  
  if (!is.null(names(group_f))) {
    common <- intersect(rownames(datExpr), names(group_f))
    if (length(common) >= 5) {
      MEs2   <- MEs[common, , drop = FALSE]
      group2 <- group_f[common]
    } else {
      warning("group_f names do not match datExpr rownames well; skipping trait export.")
      MEs2 <- NULL
    }
  } else {
    MEs2   <- MEs
    group2 <- group_f
  }
  
  if (!is.null(MEs2)) {
    trait <- data.frame(
      Control = as.numeric(group2 == "Control"),
      AsymAD  = as.numeric(group2 == "AsymAD"),
      AD      = as.numeric(group2 == "AD")
    )
    rownames(trait) <- rownames(MEs2)
    
    modTraitCor <- cor(MEs2, trait, use = "pairwise.complete.obs")
    modTraitP   <- corPvalueStudent(modTraitCor, nSamples = nrow(MEs2))
    
    write_tsv(as.data.frame(modTraitCor), file.path(out_dir, "WGCNA_ME_trait_correlation.tsv"))
    write_tsv(as.data.frame(modTraitP),   file.path(out_dir, "WGCNA_ME_trait_correlation_pvalues.tsv"))
  }
}

# -------------------- 6) Save mapping + MEs + RDS bundle ----------------------

log_step("WGCNA | save outputs")

geneModule <- data.frame(
  GeneSymbol  = colnames(datExpr),
  ModuleColor = moduleColors,
  ModuleLabel = moduleLabels,
  stringsAsFactors = FALSE
)

write_csv(geneModule, file.path(out_dir, "WGCNA_GeneModuleMap.csv"))

MEs_out <- data.frame(Sample = rownames(MEs), MEs, check.names = FALSE)
write_csv(MEs_out, file.path(out_dir, "WGCNA_ModuleEigengenes.csv"))

saveRDS(
  list(
    expr_gene     = expr_gene,
    datExpr       = datExpr,
    net           = net,
    moduleColors  = moduleColors,
    moduleLabels  = moduleLabels,
    MEs           = MEs,
    geneModule    = geneModule,
    softPower     = softPower,
    sft           = sft,
    kME_cut       = kME_cut,
    p_cut         = p_cut
  ),
  file = file.path(out_dir, "WGCNA_results.rds")
)

# -------------------- 7) Plot: dendrogram + module color bar ------------------
# ------------------------------------------------------------------------------
# Publication-quality WGCNA figure:
#   - Gene dendrogram (triangle style)
#   - Module color bar (thick, manually drawn with rect)
#   - Optional legend (top modules only)
#
# Requirements (objects must already exist):
#   - net           : result from blockwiseModules()
#   - datExpr       : samples x genes expression matrix
#   - moduleColors  : color vector for all genes in datExpr (length = ncol(datExpr))
#   - out_dir       : output directory path (e.g., "outputs")
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(WGCNA)
  library(dendextend)
})

log_step("WGCNA | publication-quality dendrogram + color bar")

if (!exists("out_dir")) out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------
# Choose the block to plot
# -------------------------------
block  <- 1
dend0  <- as.dendrogram(net$dendrograms[[block]])
bgenes <- net$blockGenes[[block]]          # indices of genes included in this block

# Module colors aligned to genes in this block
cols_leaf <- labels2colors(net$colors)[bgenes]

# Make grey less visually dominant (visual-only; does not alter results)
cols_leaf_plot <- cols_leaf
cols_leaf_plot[cols_leaf_plot == "grey"] <- "grey95"

# Remove leaf labels for readability
labels(dend0) <- rep("", length(labels(dend0)))

# Clean dendrogram styling
dend <- dend0 |>
  dendextend::set("branches_lwd", 1.2) |>
  dendextend::set("branches_col", "grey20")

dend <- dendextend::hang.dendrogram(dend, hang_height = 0)
# If you want a bit more whitespace at the bottom of triangle:
# dend <- set(dend, "hang_leaves", 0.02)

# -------------------------------
# Output settings (publication)
#   - Prefer PDF for vector graphics
#   - Also provide high-DPI PNG if needed
# -------------------------------

pdf_file <- file.path(out_dir, "FIG_WGCNA_dendrogram_moduleColors.pdf")

# Use a slightly taller page to give the color bar real space
pdf(pdf_file, width = 7.2, height = 6.8, useDingbats = FALSE)

# 3-row layout: dendrogram / colorbar / legend (separate panels)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(4.2, 1.6, 1.0))

# --- Panel 1: dendrogram
par(mar = c(1, 4.0, 2.0, 1.0), xaxs = "i", yaxs = "i")
dmax <- max(attr(dend, "height"))
plot(
  dend,
  type    = "triangle",
  leaflab = "none",
  main    = "WGCNA gene dendrogram",
  ylab    = "Height",
  ylim    = c(0, dmax * 1.15) 
)

# --- Panel 2: thick module color bar
par(mar = c(1.6, 4.0, 0.4, 1.0), xaxs = "i", yaxs = "i")
plot.new()
n <- length(cols_leaf_plot)
plot.window(xlim = c(0, n), ylim = c(0, 1))

bar_height <- 0.88
y0 <- (1 - bar_height) / 2
y1 <- y0 + bar_height

rect(
  xleft   = 0:(n - 1),
  ybottom = y0,
  xright  = 1:n,
  ytop    = y1,
  col     = cols_leaf_plot,
  border  = NA
)
axis(2, at = 0.5, labels = "Module", las = 1, tick = FALSE)
box(lwd = 0.8)

# --- Panel 3: legend (separate, no overlap)
par(mar = c(0.2, 4.0, 0.2, 1.0))
plot.new()

top_n <- 8
tab <- sort(table(cols_leaf), decreasing = TRUE)
tab <- tab[names(tab) != "grey"]

if (length(tab) > 0) {
  show <- head(tab, top_n)
  legend(
    "center",
    legend = paste0(names(show), " (n=", as.integer(show), ")"),
    fill   = names(show),
    border = NA,
    bty    = "n",
    ncol   = 2,
    cex    = 0.85
  )
}

dev.off()

cat("Saved publication-quality WGCNA figure to:", pdf_file, "\n")
cat("softPower:", softPower, "\n")
cat("nGenes:", ncol(datExpr), " nSamples:", nrow(datExpr), "\n")
cat("Module counts:\n"); print(sort(table(moduleColors), decreasing = TRUE))
cat("Grey fraction:", mean(moduleColors == "grey"), "\n")

# -------------------- 8) Heatmap: samples x hub genes -------------------------
# ------------------------------------------------------------------------------
# Uses:
#   - datExpr   : samples x genes matrix
#   - sig_genes : data.frame with GeneSymbol
#   - out_dir   : output directory
#
# Output:
#   - outputs/WGCNA_sigHubGenes_heatmap.pdf
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(pheatmap)
})

log_step("WGCNA | heatmap of significant hub genes")

heatmap_file <- file.path(out_dir, "WGCNA_sigHubGenes_heatmap.pdf")

if (!exists("sig_genes") || nrow(sig_genes) == 0) {
  warning("sig_genes is empty; skipping heatmap.")
} else {
  heat_genes <- unique(sig_genes$GeneSymbol)
  heat_genes <- intersect(heat_genes, colnames(datExpr))
  
  if (length(heat_genes) < 2) {
    warning("Not enough significant genes found in datExpr for heatmap; skipping.")
  } else {
    # Expression matrix for heatmap: genes x samples
    heat_mat <- t(datExpr[, heat_genes, drop = FALSE])
    
    # Row-wise z-score scaling for visualization
    heat_mat_scaled <- t(scale(t(heat_mat)))
    heat_mat_scaled[is.na(heat_mat_scaled)] <- 0
    
    # Optional sample annotation if group_f exists
    ann_col <- NULL
    if (exists("group_f")) {
      if (!is.null(names(group_f))) {
        grp <- group_f[colnames(heat_mat_scaled)]
      } else {
        grp <- group_f
        names(grp) <- colnames(heat_mat_scaled)
        grp <- grp[colnames(heat_mat_scaled)]
      }
      ann_col <- data.frame(Group = as.character(grp))
      rownames(ann_col) <- colnames(heat_mat_scaled)
    }
    
    # Dynamic sizing so labels remain readable
    n_genes   <- nrow(heat_mat_scaled)
    n_samples <- ncol(heat_mat_scaled)
    
    pdf(
      heatmap_file,
      width  = max(10, n_samples * 0.35 + 4),
      height = max(8,  n_genes   * 0.28 + 4),
      useDingbats = FALSE
    )
    
    pheatmap(
      heat_mat_scaled,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      annotation_col = ann_col,
      fontsize_row = ifelse(n_genes <= 40, 10, ifelse(n_genes <= 80, 8, 6)),
      fontsize_col = ifelse(n_samples <= 40, 10, ifelse(n_samples <= 80, 8, 6)),
      cellwidth  = ifelse(n_samples <= 50, 16, 10),
      cellheight = ifelse(n_genes   <= 50, 14, 10),
      border_color = "grey85",
      scale = "none",
      main = "WGCNA significant hub genes expression heatmap\n"
    )
    
    dev.off()
    cat("Saved heatmap to:", heatmap_file, "\n")
  }
}

# -------------------- 9) Horizontal barplot: module sizes ---------------------
# ------------------------------------------------------------------------------
# Uses:
#   - geneModule : data.frame with ModuleColor
#   - out_dir    : output directory
#
# Output:
#   - outputs/WGCNA_module_size_barplot.pdf
# ------------------------------------------------------------------------------

log_step("WGCNA | horizontal barplot of module sizes")

barplot_file <- file.path(out_dir, "WGCNA_module_size_barplot.pdf")

if (!exists("geneModule") || nrow(geneModule) == 0) {
  warning("geneModule is empty; skipping module barplot.")
} else {
  mod_counts <- sort(table(geneModule$ModuleColor), decreasing = FALSE)
  mod_names  <- names(mod_counts)
  
  # Make grey lighter for plotting only
  mod_fill <- mod_names
  mod_fill[mod_fill == "grey"] <- "grey80"
  
  pdf(
    barplot_file,
    width = 10,
    height = max(6, length(mod_counts) * 0.35 + 2),
    useDingbats = FALSE
  )
  
  par(mar = c(4, 10, 3, 2))
  bp <- barplot(
    mod_counts,
    horiz = TRUE,
    las = 1,
    col = mod_fill,
    border = NA,
    xlab = "Number of genes",
    main = "Number of genes in each WGCNA module\n",
    cex.names = 0.9
  )
  
  text(
    x = as.numeric(mod_counts) + max(mod_counts) * 0.01,
    y = bp,
    labels = as.numeric(mod_counts),
    pos = 4,
    cex = 0.85
  )
  
  box()
  dev.off()
  
  cat("Saved module barplot to:", barplot_file, "\n")
}

log_step("WGCNA | plots done")
