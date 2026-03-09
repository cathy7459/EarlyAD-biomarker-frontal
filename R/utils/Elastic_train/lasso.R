suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(data.table)
})

train_lasso <- function(X, y, cv_folds = 10, seed = 1) {
  set.seed(seed)
  cv <- cv.glmnet(x = as.matrix(X), y = y, family = "binomial", alpha = 1, nfolds = cv_folds)
  list(cv = cv)
}

pick_lambda <- function(cv, which = c("lambda.1se", "lambda.min")) {
  which <- match.arg(which)
  if (which == "lambda.1se") cv$lambda.1se else cv$lambda.min
}

selected_genes <- function(cv, lambda) {
  b <- as.matrix(coef(cv, s = lambda))
  nz <- which(b[,1] != 0)
  data.table::data.table(
    feature = rownames(b)[nz],
    coef = b[nz,1]
  )[feature != "(Intercept)"]
}

roc_plot_to_file <- function(prob, y, path) {
  roc_obj <- pROC::roc(y, prob, quiet = TRUE)
  png(path, width = 900, height = 700)
  plot(roc_obj, main = sprintf("ROC (AUC=%.3f)", pROC::auc(roc_obj)))
  dev.off()
  list(roc = roc_obj, auc = as.numeric(pROC::auc(roc_obj)))
}
