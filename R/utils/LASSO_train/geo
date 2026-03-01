suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

# Read only the series matrix table block (!series_matrix_table_begin ~ end)
read_series_matrix_table <- function(gzfile_path) {
  if (!file.exists(gzfile_path)) stop("Series matrix not found: ", gzfile_path)

  con <- gzfile(gzfile_path, open = "rt")
  on.exit(close(con), add = TRUE)

  line <- ""
  while (length(line) > 0) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^!series_matrix_table_begin", line)) break
  }
  if (length(line) == 0) stop("Could not find !series_matrix_table_begin in: ", gzfile_path)

  tmp <- tempfile(fileext = ".txt")
  out <- file(tmp, open = "wt")
  on.exit({ try(close(out), silent = TRUE) }, add = TRUE)

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^!series_matrix_table_end", line)) break
    writeLines(line, con = out)
  }
  close(out)

  dt <- data.table::fread(tmp)
  unlink(tmp)
  dt
}

pick_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) > 0) return(hit[1])

  # case-insensitive partial match
  cols_l <- tolower(cols)
  cand_l <- tolower(candidates)
  for (cand in cand_l) {
    idx <- which(grepl(cand, cols_l, fixed = TRUE))
    if (length(idx) > 0) return(cols[idx[1]])
  }
  NA_character_
}

get_gpl_annotation <- function(gpl_id) {
  gpl <- GEOquery::getGEO(gpl_id, AnnotGPL = TRUE)
  tab <- GEOquery::Table(gpl)
  data.table::as.data.table(tab)
}
