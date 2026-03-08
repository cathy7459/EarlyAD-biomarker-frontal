# 00_setup.R
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(glmnet)
  library(pROC)
  library(stringr)
})

source("R/utils/LASSO_test/geo.R")
source("R/utils/LASSO_test/model.R")
source("config//LASSO_test/config.R")

# output dirs
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

state <- list(cfg = CFG)
log_msg("Config loaded. prefix=", CFG$prefix)
