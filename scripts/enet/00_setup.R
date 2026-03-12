suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(glmnet)
  library(pROC)
  library(stringr)
})

#source(file.path("config", "config.R"))
#source(file.path("R", "utils_io.R"))
#source(file.path("R", "utils_geo.R"))
#source(file.path("R", "utils_model.R"))

.dir_create(CFG$out_dir)
.dir_create(file.path(CFG$out_dir, "tables"))
.dir_create(file.path(CFG$out_dir, "figures"))
.dir_create(file.path(CFG$out_dir, "logs"))
.dir_create(file.path("data", "raw"))
.dir_create(file.path("data", "processed"))

log_file <- file.path(CFG$out_dir, "logs", paste0(CFG$prefix, "_log.txt"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

if (!exists("state") || is.null(state)) state <- new.env(parent = emptyenv())
state$cfg <- CFG

log_msg("Config loaded. prefix=", CFG$prefix)
log_msg("R version detected: ", R.version.string)
