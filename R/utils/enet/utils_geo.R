# R/utils_geo.R
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
})

pick_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) > 0) return(hit[1])

  for (cand in candidates) {
    idx <- which(tolower(cols) == tolower(cand))
    if (length(idx) > 0) return(cols[idx[1]])
  }

  for (cand in candidates) {
    idx <- grep(tolower(cand), tolower(cols), fixed = TRUE)
    if (length(idx) > 0) return(cols[idx[1]])
  }

  NA_character_
}

clean_symbol <- function(x) {
  x <- gsub('"', "", as.character(x))
  x <- trimws(x)
  x <- strsplit(x, " /// |///|//|;|,")[[1]][1]
  x <- trimws(x)
  if (is.na(x) || x == "" || x == "---") return(NA_character_)
  x
}

read_series_matrix_table <- function(gzfile) {
  if (!file.exists(gzfile)) stop("Series matrix not found: ", gzfile)

  con <- gzfile(gzfile, open = "rt")
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) stop("!series_matrix_table_begin not found")
    if (grepl("^!series_matrix_table_begin", line)) break
  }

  header <- readLines(con, n = 1)
  header <- gsub('"', "", header)
  header <- trimws(strsplit(header, "\t")[[1]])

  buf <- character()
  repeat {
    x <- readLines(con, n = 5000)
    if (length(x) == 0) break
    end <- which(grepl("^!series_matrix_table_end", x))
    if (length(end) > 0) {
      buf <- c(buf, x[seq_len(end[1] - 1)])
      break
    }
    buf <- c(buf, x)
  }

  dt <- fread(text = paste(buf, collapse = "\n"), sep = "\t", header = FALSE, quote = "")
  setnames(dt, header)
  if (!("ID_REF" %in% names(dt))) stop("ID_REF column missing after read")
  dt[, ID_REF := as.character(ID_REF)]
  dt
}

parse_series_matrix_sample_meta <- function(gzfile) {
  if (!file.exists(gzfile)) stop("Series matrix not found: ", gzfile)

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

  if (is.na(id_col)) id_col <- names(gpl_tab)[1]
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
  unique(anno[!is.na(ID_REF) & ID_REF != "" & !is.na(GeneSymbol) & GeneSymbol != "", .(ID_REF, GeneSymbol)])
}

collapse_to_gene_maxvar <- function(expr_mat, map_dt, keep_genes = NULL) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat <- as.matrix(expr_mat)
  
  if (is.null(rownames(expr_mat))) {
    stop("expr_mat must have probe IDs in rownames.")
  }
  
  required_cols <- c("ProbeID", "GeneSymbol")
  if (!all(required_cols %in% names(map_dt))) {
    stop("map_dt must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  expr_probe_ids <- normalize_probe_id(rownames(expr_mat))
  rownames(expr_mat) <- expr_probe_ids
  
  map_dt <- data.table::as.data.table(map_dt)
  map_dt[, ProbeID := normalize_probe_id(ProbeID)]
  map_dt[, GeneSymbol := normalize_gene_symbol(GeneSymbol)]
  
  map_dt <- map_dt[
    !is.na(ProbeID) & ProbeID != "" &
      !is.na(GeneSymbol) & GeneSymbol != ""
  ]
  
  if (!is.null(keep_genes)) {
    keep_genes <- normalize_gene_symbol(keep_genes)
    keep_genes <- unique(keep_genes[!is.na(keep_genes) & keep_genes != ""])
    map_dt <- map_dt[GeneSymbol %in% keep_genes]
  }
  
  overlap_probe <- intersect(rownames(expr_mat), map_dt$ProbeID)
  if (length(overlap_probe) == 0) {
    stop("No overlap probes/genes after filtering")
  }
  
  expr_sub <- expr_mat[overlap_probe, , drop = FALSE]
  map_sub <- map_dt[match(overlap_probe, ProbeID), .(ProbeID, GeneSymbol)]
  
  probe_var <- apply(expr_sub, 1, var, na.rm = TRUE)
  probe_var[!is.finite(probe_var)] <- -Inf
  
  collapse_dt <- data.table::data.table(
    ProbeID = overlap_probe,
    GeneSymbol = map_sub$GeneSymbol,
    variance = probe_var
  )
  
  collapse_dt <- collapse_dt[!is.na(GeneSymbol) & GeneSymbol != ""]
  data.table::setorder(collapse_dt, GeneSymbol, -variance, ProbeID)
  best_dt <- collapse_dt[, .SD[1], by = GeneSymbol]
  
  expr_gene_mat <- expr_sub[best_dt$ProbeID, , drop = FALSE]
  rownames(expr_gene_mat) <- best_dt$GeneSymbol
  
  expr_gene_mat
}

normalize_sample_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^X+", "", x)
  x <- gsub("\\.CEL(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- gsub("\\.gz$", "", x, ignore.case = TRUE)
  x <- gsub("\\.txt$", "", x, ignore.case = TRUE)
  x <- gsub("\\.csv$", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("[\"']", "", x)
  x <- gsub("\\.", "-", x)
  toupper(x)
}
