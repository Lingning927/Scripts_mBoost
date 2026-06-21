library(mixOmics)
library(randomForest)
library(e1071)
library(Rtsne)
library(caret)
library(mBoost)
library(glmnet)
library(xgboost)
library(keras3)
library(pROC)
library(parallel)

# Example 6: negative-control simulation with random labels. It checks
# whether the diagnostics avoid reporting spurious model improvement.
get_auc <- function(true_labels, pred_probs) {
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE)
  return(as.numeric(auc(roc_obj)))
}

# Train one model type, compute train/test AUC, and estimate RP/IP statistics
# under label perturbation.
overfit_model_type <- function(X_train, Y_train, X_test, Y_test, model_type, boot_num = 100) {
  pos_class <- levels(Y_train)[2]
  if(model_type == "knn") {
    model_knn <- knn3(X_train, Y_train, k = 5)
    prob_train <- predict(model_knn, X_train, type = "prob")[, pos_class]
    prob_test <- predict(model_knn, X_test, type = "prob")[, pos_class]
    model_train <- function(X_train, Y_train) {
        model_knn <- knn3(X_train, Y_train, k = 5)
        prob_train <- predict(model_knn, X_train, type = "prob")[, pos_class]
        auc0 <- get_auc(Y_train, prob_train)
        return(auc0)
    }
  }else if(model_type == "svm") {
    model_svm <- svm(X_train, Y_train, probability = TRUE, kernel = "radial")
    pred_train_svm <- predict(model_svm, X_train, probability = TRUE)
    prob_train <- attr(pred_train_svm, "probabilities")[, pos_class]
    pred_test_svm <- predict(model_svm, X_test, probability = TRUE)
    prob_test <- attr(pred_test_svm, "probabilities")[, pos_class]
    model_train <- function(X_train, Y_train) {
        model_svm <- svm(X_train, Y_train, probability = TRUE, kernel = "radial")
        pred_train_svm <- predict(model_svm, X_train, probability = TRUE)
        prob_train <- attr(pred_train_svm, "probabilities")[, pos_class]
        auc0 <- get_auc(Y_train, prob_train)
        return(auc0)
    }
  }else if (model_type == "rf") {
    model_rf <- randomForest(x = X_train, y = Y_train, ntree = 500)
    prob_train <- predict(model_rf, X_train, type = "prob")[, pos_class]
    prob_test <- predict(model_rf, X_test, type = "prob")[, pos_class]
    model_train <- function(X_train, Y_train) {
        model_rf <- randomForest(x = X_train, y = Y_train, ntree = 500)
        prob_train <- predict(model_rf, X_train, type = "prob")[, pos_class]
        auc0 <- get_auc(Y_train, prob_train)
        return(auc0)
    }
  }else if(model_type == "xgb") {
    Y_train_num <- ifelse(Y_train == pos_class, 1, 0)
    Y_test_num <- ifelse(Y_test == pos_class, 1, 0)
    dtrain <- xgb.DMatrix(data = X_train, label = Y_train_num)
    dtest <- xgb.DMatrix(data = X_test, label = Y_test_num)
    params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2, eta = 0.01)
    model_xgb <- xgb.train(params = params, data = dtrain, nrounds = 50)
    prob_train <- predict(model_xgb, dtrain)
    prob_test <- predict(model_xgb, dtest)
    model_train <- function(X_train, Y_train) {
        Y_train_num <- ifelse(Y_train == pos_class, 1, 0)
        dtrain <- xgb.DMatrix(data = X_train, label = Y_train_num)
        model_xgb <- xgb.train(params = params, data = dtrain, nrounds = 100)
        prob_train <- predict(model_xgb, dtrain)
        auc0 <- get_auc(Y_train, prob_train)
        return(auc0)
    }
  }else if (model_type == "lasso") {
    Y_train_num <- ifelse(Y_train == pos_class, 1, 0)
    Y_test_num <- ifelse(Y_test == pos_class, 1, 0)
    cv_lasso <- cv.glmnet(X_train, Y_train_num, family = "binomial", alpha = 1)
    best_lambda <- cv_lasso$lambda.min
    prob_train <- predict(cv_lasso, newx = X_train, s = best_lambda, type = "response")[, 1]
    prob_test <- predict(cv_lasso, newx = X_test, s = best_lambda, type = "response")[, 1]
    model_train <- function(X_train, Y_train) {
        Y_train_num <- ifelse(Y_train == pos_class, 1, 0)
        cv_lasso <- cv.glmnet(X_train, Y_train_num, family = "binomial", alpha = 1)
        best_lambda <- cv_lasso$lambda.min
        prob_train <- predict(cv_lasso, newx = X_train, s = best_lambda, type = "response")[, 1]
        auc0 <- get_auc(Y_train, prob_train)
        return(auc0)
    }
  }
  auc0 <- get_auc(Y_train, prob_train)
  auc_test <- get_auc(Y_test, prob_test)
  Y_train_num <- ifelse(Y_train == pos_class, 1, 0)
  IP <- ps_probs(X_train, Y_train_num, prob_train, 20, 200)
  boot_glm <- function(boot) {
    set.seed(boot + 100)
    Y2 <- Y_train[sample(1:length(Y_train), length(Y_train))]
    acc0 <- model_train(X_train, Y2)
    return(acc0)
  }
  accs <- mclapply(seq(1, boot_num), boot_glm, mc.cores = 20)
  accs <- unlist(accs)
  P <- (sum(auc0 <= accs) + 1) / (boot_num + 1)
  return(list(acc0 = auc0, acc_test = auc_test, P = P, IP = IP))
}

args <- commandArgs(trailingOnly = TRUE)
seed <- ifelse(length(args) >= 1, as.integer(args[1]), 1)

n <- 100
p <- 20
set.seed(seed)

X_train <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
X_test <- mvrnorm(round(0.2*n), mu = rep(0, p), Sigma = diag(p))

Y <- sample(c(0, 1), nrow(X_train), replace = TRUE)
Y_train <- as.factor(Y)

Y <- sample(c(0, 1), nrow(X_test), replace = TRUE)
Y_test <- as.factor(Y)

models <- c("knn", "rf", "lasso", "xgb", "svm")
results <- matrix(0, 4, length(models))
colnames(results) <- models
rownames(results) <- c("Train_AUC", "Test_AUC", "RP", "IP")
for(i in 1:length(models)) {
    model_type <- models[i]
    tmp <- overfit_model_type(X_train, Y_train, X_test, Y_test, model_type)
    results[, i] <- as.vector(unlist(tmp))
}

saveRDS(results, paste0("simulation/simu_RP/id_", seed, ".rds"))
