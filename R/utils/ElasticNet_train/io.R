suppressPackageStartupMessages({
  library(data.table)
})

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

read_csv_dt <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  data.table::fread(path)
}
