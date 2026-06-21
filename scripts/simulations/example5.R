library(ggplot2)
library(glmnet)
library(MASS)
library(parallel)
library(randomForest)

# Example 5: nonlinear feature-selection simulation. Run this script with a
# seed and a feature-selection method, e.g. Rscript example5.R 1 RF.

args <- commandArgs(trailingOnly = TRUE)
seed <- ifelse(length(args) >= 1, as.integer(args[1]), 1)
method <- ifelse(length(args) >= 2, args[2], "RF")

true_model <- function(X) {
    # Only the first two variables determine the label; all others are noise.
    x1 <- X[, 1]
    x2 <- X[, 2]
    Y <- as.numeric((x1^2 + x2^2) < 8 / pi)
    return(Y)
}

gene_X <- function(n, p){
    # Keep the first two features informative and fill the rest with noise.
    x1 <- runif(n, -2, 2)
    x2 <- runif(n, -2, 2)
    X <- mvrnorm(n, mu = rep(0, p), Sigma = 8*diag(p))
    X[, 1] <- x1
    X[, 2] <- x2
    return(X)
}

# Apply one of the feature-selection methods used in the manuscript.
my_fs_simu <- function(X, Y, method = "RF") {
  if (method == "RF") {
    rf_model <- randomForest(X, as.factor(Y), importance=TRUE)
    top_var <- order(rf_model$importance[, "MeanDecreaseAccuracy"], decreasing = TRUE)[1:10]
  }else if(method == "Lasso") {
    fit <- cv.glmnet(X, Y, family = "binomial")
    beta_hat <- fit$glmnet.fit$beta[, fit$index[1]]
    top_var <- which(beta_hat != 0)
  }else if(method == "Wilcoxon") {
    significant_features <- c()
    scores <- c()
    for (i in 1:ncol(X)) {
      group1 <- X[Y == 1, i]
      group2 <- X[Y == 0, i]
      test_result <- wilcox.test(group1, group2)
      if (is.na(test_result$p.value)) {
          scores <- c(scores, 1)
          next
      }else if (test_result$p.value < 0.05) {
        significant_features <- c(significant_features, i)
        scores <- c(scores, test_result$p.value)
      }else {
        scores <- c(scores, test_result$p.value)
      }
    }
    top_var <- significant_features
    if(length(significant_features) < 5) {
        top_var <- order(scores)[1:5]
    }
  }else if(method == "KS") {
    significant_features <- c()
    for (i in 1:ncol(X)) {
      test_result <- ks.test(X[, i] ~ Y)
      if (is.na(test_result$p.value)) {
          next
      }else if (test_result$p.value < 0.05) {
        significant_features <- c(significant_features, i)
      }
    }
    top_var <- significant_features
  }else if(method == "T_test") {
    significant_features <- c()
    scores <- c()
    for (i in 1:ncol(X)) {
      test_result <- t.test(X[, i] ~ Y)
      if (is.na(test_result$p.value)) {
          scores <- c(scores, 1)
          next
      }else if (test_result$p.value < 0.05) {
        significant_features <- c(significant_features, i)
        scores <- c(scores, test_result$p.value)
      }else {
        scores <- c(scores, test_result$p.value)
      }
    }
    top_var <- significant_features
    if(length(significant_features) < 5) {
        top_var <- order(scores)[1:5]
    }
    top_var <- significant_features
  }
  return(top_var)
}

# Feature-validity p-value based on overlap stability across random splits.
pv_simu <- function(predictors, response, method, reps = 10, BootNum = 100) {
    sub_part <- function(seed, id_pre, predictors, response, method) {
        if(seed > 0) {
            set.seed(seed)
            m <- sum(response == 1)
            id <- sample(seq(1, length(response)), m)
            response2 <- rep(0, length(response))
            response2[id] <- 1
        }else {
            response2 <- response
        }
        f1 <- my_fs_simu(predictors[id_pre, ], response2[id_pre], method)
        f2 <- my_fs_simu(predictors[-id_pre, ], response2[-id_pre], method)
        return(length(intersect(f1, f2)))
    }
    p_vals <- rep(0, reps)
    for (i in 1:reps) {
        set.seed(i+1000)
        id_pre <- sample(seq(1, length(response)), round(length(response) / 2))
        fs <- mclapply(seq(1, BootNum), FUN = sub_part, id_pre, predictors,
                response, method, mc.cores = 10)
        fs <- do.call(rbind, fs)
        f0 <- sub_part(0, id_pre, predictors, response, method)
        p_vals[i] <- (sum(fs > f0) + 1) / (BootNum + 1)
    }
    return(p_vals)
}

# Feature-redundancy p-value for a Random Forest trained on selected features.
rf_rp_compute <- function(X, Y, boot_num = 100) {
  sub_rf <- function(seed, X, Y) {
    if(seed > 0) {
      set.seed(seed)
      Y2 <- Y[sample(seq(1, n), n)]
    }else {
      Y2 <- Y
    }
    rf_model <- randomForest(X, as.factor(Y2), maxnodes = 5)
    acc <- sum(Y2 == predict(rf_model, X)) / length(Y2)
    return(acc)
  }
  accs <- mclapply(seq(1, boot_num), sub_rf, X, Y, mc.cores = 10)
  accs <- unlist(accs)
  acc0 <- sub_rf(0, X, Y)
  Pr <- (sum(acc0 <= accs) + 1) / (boot_num + 1)
  return(Pr)
}

n <- 200
p_list <- seq(100, 500, 100)

res <- matrix(0, length(p_list), 4)
i <- 1
for (p in p_list) {
    set.seed(seed)
    X <- gene_X(n, p)
    Y <- true_model(X)

    X_test <- gene_X(round(0.2*n), p)
    Y_test <- true_model(X_test)

    if (method == "All") {
        rf_model <- randomForest(X, as.factor(Y), importance=TRUE, maxnodes = 5)
        acc_train <- sum(Y == predict(rf_model, X)) / length(Y)
        acc_test <- sum(Y_test == predict(rf_model, X_test)) / length(Y_test)
        Pr <- rf_rp_compute(X, Y)
        res[i, ] <- c(acc_train, acc_test, Pr, 0)
    }else {
        top_var <- my_fs_simu(X, Y, method)
        X2 <- X[, top_var]
        rf_model <- randomForest(X2, as.factor(Y), importance=TRUE, maxnodes = 5)
        acc_train <- sum(Y == predict(rf_model, X[, top_var])) / length(Y)
        acc_test <- sum(Y_test == predict(rf_model, X_test[, top_var])) / length(Y_test)
        Pr <- rf_rp_compute(X2, Y)
        Pv <- mean(pv_simu(X, Y, method))
        res[i, ] <- c(acc_train, acc_test, Pr, Pv)
    }
    i <- i + 1
}
rownames(res) <- as.character(p_list)
colnames(res) <- c("acc_train", "acc_test", "Pr", "Pv")
