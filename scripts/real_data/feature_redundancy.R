require(SIAMCAT)
require(parallel)
require(dplyr)

# Train a SIAMCAT classifier and return its internal cross-validation AUC.
# When boot_seed > 0, labels are permuted to form the null distribution for pr.
compute_auc <- function(feat, meta, boot_seed = 0, method = "lasso") {
  if (boot_seed > 0) {
    set.seed(100 + boot_seed)
    n <- dim(meta)[1]
    meta$Group <- meta$Group[sample(1:n, n)]
  }
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
      method = "lasso",
      feature.type = "normalized",
      verbose = 0
    )
  siamcat.train <- make.predictions(siamcat.train, verbose = 0)
  siamcat.train <-  evaluate.predictions(siamcat.train, verbose = 0)

  auc <- as.numeric(eval_data(siamcat.train)$auroc)
  return(auc)
}

# Feature-redundancy p-value from the observed AUC versus permuted-label AUCs.
compute_aucs <- function(feat, meta, boot_num = 200) {
  auc0 <- compute_auc(feat, meta)
  aucs <- mclapply(1:boot_num, FUN = compute_auc, feat = feat,
    meta = meta, mc.cores = 10)
  aucs <- unlist(aucs)
  P <- (sum(auc0 < aucs) + 1) / (boot_num + 1)
  return(list(auc0 = auc0, aucs = aucs, P = P))
}

# For each cohort, compute pr and cross-cohort AUCs against other cohorts of
# the same disease/data type.
Pr_one_cohort <- function(feat_list, meta_list, method = "lasso", boot_num = 200) {
  project.list <- names(feat_list)
  oto_res <- list()
  for (project in project.list) {
    feat <- feat_list[[project]]
    meta <- meta_list[[project]]
    all_feacture <- colnames(feat)
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
    aucs <- mclapply(1:boot_num, FUN = compute_auc, feat = feat,
      meta = meta, mc.cores = 10)
    aucs <- unlist(aucs)
    P <- (sum(auc_train < aucs) + 1) / (boot_num + 1)
    temp_list <- list()
    temp_list[["auc_train"]] <- auc_train
    temp_list[["Pr"]] <- P
    for (project2 in project.list) {
      if (project == project2) {
        next
      }
      feat_test <- feat_list[[project2]]
      meta_test <- meta_list[[project2]]
      lack_feactures <- setdiff(all_feacture, colnames(feat_test))
      tmp <- matrix(0, nrow = dim(feat_test)[1],
                    ncol = length(lack_feactures))
      colnames(tmp) <- lack_feactures
      tmp <- cbind(feat_test, tmp)
      feat_test <- tmp[, all_feacture]
      siamcat.test <- siamcat(feat = t(feat_test), meta=as.data.frame(meta_test),
                            label='Group', case='Case', verbose = 0)
      
      siamcat.test <- normalize.features(siamcat.test,
                                          norm.param=norm_params(siamcat.train),
                                          feature.type = 'original',
                                          verbose = 0)
      siamcat.test <- make.predictions(
        siamcat = siamcat.train,
        siamcat.holdout = siamcat.test,
        normalize.holdout = TRUE, verbose = 0)
      siamcat.test <- evaluate.predictions(siamcat.test, verbose = 0)
      auc_test <- as.numeric(eval_data(siamcat.test)$auroc)
      temp_list[[project2]] <- auc_test
    }
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
res_pr <- Pr_one_cohort(feat_list, meta_list, method)
