library(SIAMCAT)
library(glmnet)
library(randomForest)
library(pROC)
library(parallel)
library(kSamples)


my_fs <- function(response, predictors, method = "Wilcoxon", fs_num = 20) {
    significant_features <- c()
    ranked_feacture <- c()
    if (method == "T_test") {
    scores <- c()
    for (i in 1:ncol(predictors)) {
        t_result <- t.test(predictors[, i] ~ response)
        if (is.na(t_result$p.value)) {
        next
        }else if (t_result$p.value < 0.05) {
        significant_features <- c(significant_features, i)
        scores <- c(scores, t_result$p.value)
        }
    }
    ranked_feacture <- colnames(predictors)[significant_features[order(scores)]]
    }

    if (method == "AD") {
    
    scores <- c()
    for (i in 1:ncol(predictors)) {
        try({
            t_result <- ad.test(predictors[, i] ~ response)
        if (is.na(t_result$ad[1, 3])) {
        next
        }else {
        significant_features <- c(significant_features, i)
        scores <- c(scores, t_result$ad[1, 3])
        }
        }, silent = TRUE)
    }
    significant_features <- significant_features[order(scores)[1:fs_num]]
    ranked_feacture <- colnames(predictors)[significant_features]
    }


    if (method == "KS") {
    significant_features <- c()
    scores <- c()
    for (i in 1:ncol(predictors)) {
        ks_result <- ks.test(predictors[, i] ~ response)
        if (is.na(ks_result$p.value)) {
        next
        }else {
        significant_features <- c(significant_features, i)
        scores <- c(scores, ks_result$p.value)
        }
    }
    significant_features <- significant_features[order(scores)[1:fs_num]]
        ranked_feacture <- colnames(predictors)[significant_features]
    }
    if (method == "Wilcoxon") {
        significant_features <- c()
        scores <- c()
        for (i in 1:ncol(predictors)) {
            group1 <- predictors[response == levels(response)[1], i]
            group2 <- predictors[response == levels(response)[2], i]
            test_result <- wilcox.test(group1, group2)
            if (is.na(test_result$p.value)) {
                next
            }else {
            significant_features <- c(significant_features, i)
            scores <- c(scores, test_result$p.value)
            }
        }
        significant_features <- significant_features[order(scores)[1:fs_num]]
        ranked_feacture <- colnames(predictors)[significant_features]
    }


    if (method == "Lasso") {
    predictors <- as.matrix(predictors)
    cv_fit <- cv.glmnet(predictors, response, family = "binomial", nfolds = 5)
    fit <- cv_fit$glmnet.fit
    beta <- fit$beta[, cv_fit$index[1]]
    id <- which(beta != 0)
    significant_features <- names(fit$beta[, cv_fit$index[1]])[id]
    ranked_feacture <- names(beta[id])[order(abs(beta[id]), decreasing = TRUE)]
    }
    if (method == "Dacc") {
        rf_model <- randomForest(x = predictors,
            y = response, importance=TRUE, ntree = 500)
        imp <- importance(rf_model)
        importance_dac <- imp[, "MeanDecreaseAccuracy"]
        significant_features <- names(importance_dac)[order(importance_dac, decreasing = TRUE)][1:fs_num]
        ranked_feacture <- significant_features
    }
    return(c(ranked_feacture))
}

feature_validaty <- function(predictors, response, method, reps = 10, BootNum = 100) {
    sub_part <- function(seed, id_pre, predictors, response, method) {
        if(seed > 0) {
            set.seed(seed)
            m <- sum(response == levels(response)[1])
            id <- sample(seq(1, length(response)), m)
            response2 <- rep(levels(response)[2], length(response))
            response2[id] <- levels(response)[1]
            response2 <- as.factor(response2)
        }else {
            response2 <- response
        }
        f1 <- my_fs(response2[id_pre], predictors[id_pre, ], method)
        f2 <- my_fs(response2[-id_pre], predictors[-id_pre, ], method)
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

pv_one_cohort <- function(feat_list, meta_list, method, reps = 10, BootNum = 100) {
  res <- NULL
  for (project in names(feat_list)) {
    feat <- feat_list[[project]]
    meta <- meta_list[[project]]
    response <- as.factor(meta$Group)
    predictors <- feat
    pvals <- feature_validaty(predictors, response, method, reps, BootNum)
    res[[project]] <- mean(pvals)
  }
  return(res)
}

disease <- "CD"
data_type <- "Amplicon_genus"
method <- "Wilcoxon"

load(paste0("data/feat_meta/", disease, "_", data_type, ".RData"))
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list

res_pv <- pv_one_cohort(feat_list, meta_list, method)
