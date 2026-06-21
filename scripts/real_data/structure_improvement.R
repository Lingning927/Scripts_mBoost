source("methods/NewMac.R")
require(parallel)
require(SIAMCAT)
require(dplyr)

# Structure-improvement p-value for SIAMCAT models, using model-predicted
# probabilities to generate bootstrap labels.
ps_based_siamcat <- function(X, siamcat.train, BootNumber = 1000) {
  predicted_probs <- eval_data(siamcat.train)$roc$predictor
  Y <- eval_data(siamcat.train)$roc$response
  Y[Y == -1] <- 0
  boot_tb <- function(boot, X, Y, predicted_probs) {
    set.seed(boot + 100)
    Y2 <- rbinom(length(Y), 1, predicted_probs)
    set.seed(boot + 101)
    Y22 <- rbinom(length(Y), 1, predicted_probs)
    mac1 <- MAC(X, Y, X, Y2)
    mac2 <- MAC(X, Y22, X, Y2)
    return(c(mac1, mac2))
  }
  Tb <- mclapply(seq(1, BootNumber), FUN = boot_tb, X,
    Y, predicted_probs, mc.cores = 10)
  Tb <- do.call(rbind, Tb)
  T0 <- mean(Tb[, 1])
  null_t <- Tb[, 2]
  P <- (sum(null_t >= T0) + 1) / (BootNumber + 1)
  return(P)
}

# Train one SIAMCAT model per cohort and summarize AUC with the ps diagnostic.
my_siamcat_with_ps <- function(feat_list, meta_list, method = "lasso", boot_num = 200) {
  project.list <- names(feat_list)
  oto_res <- list()
  for (project in project.list) {
    feat <- feat_list[[project]]
    meta <- meta_list[[project]]
    siamcat.train <- siamcat(feat = t(feat), meta=as.data.frame(meta),
                            label='Group', case='Case', verbose = 0)
    siamcat.train <- filter.features(
      siamcat.train,
      filter.method = 'abundance',
      cutoff = 0.001,
      rm.unmapped = TRUE,
      verbose=0
    )
    siamcat.train <- normalize.features(
      siamcat.train,
      norm.method = "log.std",
      norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1),
      verbose = 0
    )
    siamcat.train <-  create.data.split(
      siamcat.train,
      num.folds = 5,
      num.resample = 3,
      verbose = 0
    )
    siamcat.train<- train.model(
      siamcat.train,
      method = method,
      feature.type = "normalized",
      verbose = 0
    )
    siamcat.train <- make.predictions(siamcat.train, verbose = 0)
    siamcat.train <-  evaluate.predictions(siamcat.train, verbose = 0)
    auc_train <- as.numeric(eval_data(siamcat.train)$auroc)
    P <- ps_based_siamcat(feat, siamcat.train)
    temp_list <- list()
    temp_list[["auc_train"]] <- auc_train
    temp_list[["Ps"]] <- P
    oto_res[[project]] <- temp_list
  }
  return(oto_res)
}

disease <- "CD"
data_type <- "Amplicon_genus"
method <- "lasso"

load(paste0("data/feat_meta/", disease, "_", data_type, ".RData"))
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list
res_ps <- my_siamcat_with_ps(feat_list, meta_list, method)
