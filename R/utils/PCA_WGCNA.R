# PCA_WGCNA.R ----
# Small helpers to keep scripts readable and reproducible.

log_step <- function(title) {
  cat("\n", paste0("## ", title), "\n", sep = "")
}

require_objects <- function(objs) {
  missing <- objs[!vapply(objs, exists, logical(1), inherits = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required object(s): ", paste(missing, collapse = ", "),
         "\nProvide them in your upstream preprocessing script, then re-run.")
  }
}

as_numeric_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

safe_scale <- function(X) {
  # Scale columns; keep matrix shape
  scale(X, center = TRUE, scale = TRUE)
}

write_lines <- function(x, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  writeLines(x, con = file, useBytes = TRUE)
}

write_csv <- function(df, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file = file, row.names = FALSE)
}

write_tsv <- function(df, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, col.names = NA)
}
