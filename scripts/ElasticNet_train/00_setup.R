suppressPackageStartupMessages({
  library(data.table)
})

source("config/ElasticNet_train.R")
source("R/utils/ElasticNet_train/io.R")
source("R/utils/ElasticNet_train/geo.R")
source("R/utils/ElasticNet_train/preprocess.R")
source("R/utils/ElasticNet_train/deg.R")
source("R/utils/ElasticNet_train/lasso.R")

ensure_dir(DIR_DATA_PROCESSED)
ensure_dir(DIR_OUT_TABLES)
ensure_dir(DIR_OUT_FIGURES)
ensure_dir(DIR_OUT_LOGS)

log_msg("Setup complete.")
