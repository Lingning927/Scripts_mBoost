library(Rtsne)
library(parallel)

coverage_ratio <- function(x, y) {
  dist_x <- as.matrix(stats::dist(x))
  diag(dist_x) <- rep(100000, dim(x)[1])
  k <- round(0.02 * length(x)) + 1
  nearest_dist_x <- apply(dist_x, 1, function(xi) {sort(xi)[k]})
  nearest_dist_xy <- apply(x, 1, function(x_row) {
    min(sqrt(colSums((t(y) - x_row)^2)))
  })
  return(sum(nearest_dist_xy < nearest_dist_x) / length(nearest_dist_xy))
}

CR_one_cohort <- function(feat_list, meta_list) {
  project.list <- names(feat_list)
  res <- list()
  for (project in project.list) {
    temp_list <- NULL
    X_train <- feat_list[[project]]
    feature <- colnames(X_train)
    for (project2 in project.list) {
      if (project == project2) {
        next
      }
      X_test <- feat_list[[project2]]
      lack_feactures <- setdiff(feature, colnames(X_test))
      tmp <- matrix(0, nrow = dim(X_test)[1],
                    ncol = length(lack_feactures))
      colnames(tmp) <- lack_feactures
      X_test <- cbind(X_test, tmp)
      X_test <- X_test[, feature]
      X_all <- as.matrix(rbind(X_train, X_test))
      tsne <- Rtsne(X_all, check_duplicates = FALSE)
      x <- tsne$Y[1:nrow(X_train), ]
      x2 <- tsne$Y[(nrow(X_train)+1):nrow(X_all), ]
      cr <- coverage_ratio(x2, x)
      temp_list[[project2]] <- cr
    }
    res[[project]] <- temp_list
  }
  return(res)
}

disease <- "CD"
data_type <- "Amplicon_genus"

load(paste0("data/feat_meta/", disease, "_", data_type, ".RData"))
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list

res_cr <- CR_one_cohort(feat_list, meta_list)
