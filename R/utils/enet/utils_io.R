# R/utils_io.R
suppressPackageStartupMessages({
  library(data.table)
})

.dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  if (exists("log_file", inherits = TRUE) && is.character(log_file) && length(log_file) == 1) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

read_csv_dt <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  data.table::fread(path)
}

require_objects <- function(objs, env = parent.frame()) {
  missing <- objs[!vapply(objs, exists, logical(1), envir = env)]
  if (length(missing) > 0) stop("Missing required objects: ", paste(missing, collapse = ", "))
}

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
}
