## Auto-extracted & cleaned utilities from original script

.dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
.dir_create(CFG$out_dir)
.dir_create(file.path(CFG$out_dir, "tables"))
.dir_create(file.path(CFG$out_dir, "figures"))
.dir_create(file.path(CFG$out_dir, "logs"))

log_file <- file.path(CFG$out_dir, "logs", paste0(CFG$prefix, "_log.txt"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

pick_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) > 0) return(hit[1])
  # case-insensitive exact match
  for (cand in candidates) {
    idx <- which(tolower(cols) == tolower(cand))
    if (length(idx) > 0) return(cols[idx[1]])
  }
  # fuzzy contains match
  for (cand in candidates) {
    idx <- grep(tolower(cand), tolower(cols), fixed = TRUE)
    if (length(idx) > 0) return(cols[idx[1]])
  }
  NA_character_
}

clean_symbol <- function(x) {
  x <- gsub('"', "", as.character(x))
  x <- trimws(x)
  # common separators in GPL tables
  x <- strsplit(x, " /// |///|//|;|,")[[1]][1]
  x <- trimws(x)
  if (is.na(x) || x == "" || x == "---") return(NA_character_)
  x
}

read_series_matrix_table <- function(gzfile) {
  con <- gzfile(gzfile, open = "rt")
  on.exit(close(con), add = TRUE)

  # Move to begin
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) stop("!series_matrix_table_begin not found")
    if (grepl("^!series_matrix_table_begin", line)) break
  }

  # Header line
  header <- readLines(con, n = 1)
  header <- gsub('"', "", header)
  header <- trimws(strsplit(header, "\t")[[1]])

  # Read until end
  buf <- character()
  repeat {
    x <- readLines(con, n = 5000)
    if (length(x) == 0) break
    end <- which(grepl("^!series_matrix_table_end", x))
    if (length(end) > 0) {
      buf <- c(buf, x[seq_len(end[1] - 1)])
      break
    } else {
      buf <- c(buf, x)
    }
  }

  dt <- fread(text = paste(buf, collapse = "\n"), sep = "\t", header = FALSE, quote = "")
  setnames(dt, header)
  if (!("ID_REF" %in% names(dt))) stop("ID_REF column missing after read")
  dt[, ID_REF := as.character(ID_REF)]
  dt
}

parse_series_matrix_sample_meta <- function(gzfile) {
  con <- gzfile(gzfile, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con)

  sample_lines <- lines[grepl("^!Sample_", lines)]
  if (length(sample_lines) == 0) stop("No !Sample_ lines found in series matrix")

  split_line <- function(line) {
    parts <- strsplit(line, "\t")[[1]]
    key <- sub("^!Sample_", "", parts[1])
    vals <- parts[-1]
    list(key = key, vals = vals)
  }

  parsed <- lapply(sample_lines, split_line)
  keys <- vapply(parsed, `[[`, character(1), "key")

  idx <- which(keys == "geo_accession")
  if (length(idx) == 0) stop("No !Sample_geo_accession line found")
  gsm <- parsed[[idx[1]]]$vals

  meta <- data.table(GSM = gsm)
  for (item in parsed) {
    if (length(item$vals) == length(gsm)) meta[[item$key]] <- item$vals
  }
  meta
}

get_probe_gene_mapping <- function(platform_id) {
  gpl_obj <- getGEO(platform_id, AnnotGPL = TRUE)
  gpl_tab <- as.data.table(Table(gpl_obj))

  # Some GPL tables are weird (header is first row). Detect and fix.
  row1 <- as.character(unlist(gpl_tab[1, ], use.names = FALSE))
  row1_lc <- tolower(row1)
  looks_like_header <- any(grepl("^id$|id_ref|probe", row1_lc)) && any(grepl("symbol|gene", row1_lc))
  if (looks_like_header) {
    setnames(gpl_tab, make.names(row1, unique = TRUE))
    gpl_tab <- gpl_tab[-1]
  }

  id_candidates  <- c("ID", "ID_REF", "Probe", "PROBE_ID", "probe_id")
  sym_candidates <- c("Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "GeneSymbol", "gene_assignment", "gene symbol", "gene_symbol")

  id_col  <- pick_col(names(gpl_tab), id_candidates)
  sym_col <- pick_col(names(gpl_tab), sym_candidates)

  if (is.na(id_col))  id_col <- names(gpl_tab)[1]
  if (is.na(sym_col)) {
    sym_hits <- grep("symbol|gene", tolower(names(gpl_tab)), value = TRUE)
    if (length(sym_hits) == 0) stop("Cannot find GeneSymbol-like column in GPL")
    sym_col <- sym_hits[1]
  }

  anno <- gpl_tab[, .(
    ID_REF = as.character(get(id_col)),
    GeneSymbol_raw = as.character(get(sym_col))
  )]
  anno[, GeneSymbol := vapply(GeneSymbol_raw, clean_symbol, character(1))]
  anno <- unique(anno[!is.na(ID_REF) & ID_REF != "" & !is.na(GeneSymbol) & GeneSymbol != "", .(ID_REF, GeneSymbol)])
  anno
}

collapse_to_gene_maxvar <- function(expr_probe_mat, anno_dt, keep_genes = NULL) {
  # expr_probe_mat: probe x sample numeric matrix
  stopifnot(!is.null(rownames(expr_probe_mat)))

  anno2 <- anno_dt[ID_REF %in% rownames(expr_probe_mat)]
  if (!is.null(keep_genes)) {
    anno2 <- anno2[GeneSymbol %in% keep_genes]
  }
  if (nrow(anno2) == 0) stop("No overlap probes/genes after filtering")

  sub_expr <- expr_probe_mat[anno2$ID_REF, , drop = FALSE]
  vars <- apply(sub_expr, 1, var, na.rm = TRUE)

  anno2 <- copy(anno2)
  anno2[, probe_var := vars[match(ID_REF, rownames(sub_expr))]]
  setorder(anno2, GeneSymbol, -probe_var)
  best <- anno2[, .SD[1], by = GeneSymbol]

  expr_gene <- expr_probe_mat[best$ID_REF, , drop = FALSE]
  rownames(expr_gene) <- best$GeneSymbol
  expr_gene <- expr_gene[sort(rownames(expr_gene)), , drop = FALSE]
  list(expr_gene = expr_gene, best_map = best)
}

standardize_by_train <- function(X_train_gxs, X_new_gxs) {
  common <- intersect(rownames(X_train_gxs), rownames(X_new_gxs))
  if (length(common) < 10) stop("Too few common genes between train and new data")

  X_train <- X_train_gxs[common, , drop = FALSE]
  X_new   <- X_new_gxs[common, , drop = FALSE]

  mu <- rowMeans(X_train, na.rm = TRUE)
  sdv <- apply(X_train, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1

  Z_train <- sweep(sweep(X_train, 1, mu, "-"), 1, sdv, "/")
  Z_new   <- sweep(sweep(X_new,   1, mu, "-"), 1, sdv, "/")

  list(Z_train = Z_train, Z_new = Z_new, mu = mu, sd = sdv, common_genes = common)
}

make_labels_gse33000 <- function(meta_dt) {
  fields <- setdiff(names(meta_dt), "GSM")
  row_text <- apply(meta_dt[, ..fields], 1, function(x) paste(na.omit(as.character(x)), collapse = " | "))

  is_hd  <- grepl("(huntington|\\bHD\\b)", row_text, ignore.case = TRUE, perl = TRUE)
  is_ad  <- grepl("(alzheimer|\\bAD\\b|load|dementia)", row_text, ignore.case = TRUE, perl = TRUE)
  is_ctl <- grepl("(control|normal|healthy|non\\s*-?dement)", row_text, ignore.case = TRUE, perl = TRUE)

  y <- rep(NA_character_, nrow(meta_dt))
  y[is_ad  & !is_ctl & !is_hd] <- "AD"
  y[is_ctl & !is_ad  & !is_hd] <- "Control"
  y[is_hd] <- "HD"

  data.table(GSM = meta_dt$GSM, y = y, row_text = row_text)
}

map_expr_to_gene_maxvar <- function(expr_feature_x_sample, platform_id, keep_genes = NULL) {
  anno <- get_probe_gene_mapping(platform_id)

  # GEOquery getGEO expr rows are feature IDs (probe IDs)
  expr_probe <- expr_feature_x_sample
  # ensure IDs are character
  rownames(expr_probe) <- as.character(rownames(expr_probe))

  collapsed <- collapse_to_gene_maxvar(expr_probe, anno, keep_genes)
  collapsed
}

