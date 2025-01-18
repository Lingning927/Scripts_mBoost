library(MASS)
library(glmnet)
library(parallel)
library(ggplot2)
library(ggsci)
library(patchwork)
source("methods/simulation_methods.R")

#The script for conducting one experiment is only provided here
#100 repeated experiments were conducted in the article

###Example logistic
n <- 200
p <- 100
real_p <- 20

set.seed(123)

gene_data <- function(n, p, real_p) {
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))  # 生成特征矩阵
  true_coef <- rnorm(p)  # 生成真实系数
  zero_index <- sample(1:p, p - real_p)
  true_coef[zero_index] <- 0
  prob <- 1 / (1 + exp(-X %*% true_coef - rnorm(n, 0, 4)))  # 计算概率
  Y <- rbinom(n, size = 1, prob = prob)
  X_test <- mvrnorm(round(0.2 * n), mu = rep(0, p), Sigma = diag(p))
  prob_test <- 1 / (1 + exp(-X_test %*% true_coef))
  Y_test <- rbinom(round(0.2 * n), size = 1, prob = prob_test)
  return(list(X = X, Y = Y, X_test = X_test, Y_test = Y_test, true_coef = true_coef))
}

lasso_fs <- function(Y, X) {
  fit <- cv.glmnet(X, Y, family = "binomial")
  beta_hat <- fit$glmnet.fit$beta[, fit$index[1]]
  est_index <- which(beta_hat != 0)
  return(est_index)
}

pv_lasso <- function(predictors, response, reps = 10, BootNum = 100) {
    sub_part <- function(seed, id_pre, predictors, response) {
        if(seed > 0) {
            set.seed(seed)
            m <- sum(response == 1)
            id <- sample(seq(1, length(response)), m)
            response2 <- rep(0, length(response))
            response2[id] <- 1
        }else {
            response2 <- response
        }
        f1 <- lasso_fs(response2[id_pre], predictors[id_pre, ])
        f2 <- lasso_fs(response2[-id_pre], predictors[-id_pre, ])
        return(length(intersect(f1, f2)))
    }
    p_vals <- rep(0, reps)
    for (i in 1:reps) {
        set.seed(i+1000)
        id_pre <- sample(seq(1, length(response)), round(length(response) / 2))
        fs <- mclapply(seq(1, BootNum), FUN = sub_part, id_pre, predictors,
                response, mc.cores = 10)
        fs <- do.call(rbind, fs)
        f0 <- sub_part(0, id_pre, predictors, response)
        p_vals[i] <- (sum(fs > f0) + 1) / (BootNum + 1)
    }
    return(p_vals)
}


my_regression <- function(data, boot_num = 200) {
  X <- data$X
  Y <- data$Y
  my_data <- data.frame(cbind(X, Y))
  fit <- glm(Y~., data = my_data, family = "binomial")
  predicted_Y <- fit$fitted.values > 0.5
  acc0 <- sum(predicted_Y == Y) / length(Y)
  test_data <- data.frame(data$X_test)
  colnames(test_data) <- colnames(my_data)[1:dim(data$X_test)[2]]
  predicted_Y_test <- predict(fit, newdata = test_data, type = "response") > 0.5
  acc_test <- sum(predicted_Y_test == data$Y_test) / length(data$Y_test)
  boot_glm <- function(boot) {
    set.seed(boot + 100)
    Y2 <- Y[sample(1:length(Y), length(Y))]
    data <- data.frame(cbind(X, Y2))
    fit <- glm(Y2 ~ ., data = data, family = "binomial")
    predicted_Y <- predict(fit, type = "response") > 0.5
    acc <- sum(predicted_Y == Y2) / length(Y2)
    return(acc)
  }
  accs <- mclapply(seq(1, boot_num), boot_glm, mc.cores = 10)
  accs <- unlist(accs)
  P_glm <- (sum(acc0 <= accs) + 1) / (boot_num + 1)
  glm_res <- list(acc = acc0, acc_test = acc_test, Pr = P_glm)

  fit <- cv.glmnet(X, Y, family = "binomial")
  beta_hat <- fit$glmnet.fit$beta[, fit$index[1]]
  true_coef <- data$true_coef
  true_index <- which(true_coef != 0)
  est_index <- which(beta_hat != 0)
  pv <- mean(pv_lasso(X, Y))
  predicted_P <- 1 / (1 + exp(- X %*% beta_hat - fit$glmnet.fit$a0[fit$index[1]]))
  predicted_Y <- as.numeric(predicted_P >= 0.5)
  acc0 <- sum(predicted_Y == Y) / length(Y)
  predicted_Y_test <- (1 / (1 + exp(- data$X_test %*% beta_hat - fit$glmnet.fit$a0[fit$index[1]]))) > 0.5
  acc_test <- sum(predicted_Y_test == data$Y_test) / length(data$Y_test)
  boot_glmnet <- function(boot) {
    set.seed(boot + 100)
    Y2 <- Y[sample(1:length(Y), length(Y))]
    fit <- cv.glmnet(X, Y2, family = "binomial")
    beta_hat <- fit$glmnet.fit$beta[, fit$index[1]]
    predicted_P <- 1 / (1 + exp(- X %*% beta_hat - fit$glmnet.fit$a0[fit$index[1]]))
    predicted_Y <- as.numeric(predicted_P >= 0.5)
    acc <- sum(predicted_Y == Y2) / length(Y2)
    return(acc)
  }
  accs <- mclapply(seq(1, boot_num), boot_glmnet, mc.cores = 10)
  accs <- unlist(accs)
  P_glmnet <- (sum(acc0 <= accs) + 1) / (boot_num + 1)
  glmnet_res <- list(acc = acc0, acc_test = acc_test, Pr = P_glmnet,
    true_index = true_index, est_index = est_index, Pv = pv)
  return(list(glm_res = glm_res, glmnet_res = glmnet_res))
}


set.seed(100)
data <- gene_data(100, 50, real_p)

df_glm <- data.frame(value = numeric(), p = numeric(), type = character())
df_glmnet <- data.frame(value = numeric(), p = numeric(), type = character())

for (p in seq(40, 120, 10)) {
    data <- gene_data(n, p, real_p)
    res <- my_regression(data)
    glm_res <- res$glm_res
    glmnet_res <- res$glmnet_res
    df_glm <- rbind(df_glm, data.frame(value = glm_res$acc, p = p, type = "Train Accuracy"))
    df_glm <- rbind(df_glm, data.frame(value = glm_res$acc_test, p = p, type = "Test Accuracy"))
    df_glm <- rbind(df_glm, data.frame(value = glm_res$Pr, p = p, type = "Pr"))
    df_glmnet <- rbind(df_glmnet, data.frame(value = glmnet_res$acc, p = p, type = "Train Accuracy"))
    df_glmnet <- rbind(df_glmnet, data.frame(value = glmnet_res$acc_test, p = p, type = "Test Accuracy"))
    df_glmnet <- rbind(df_glmnet, data.frame(value = glmnet_res$Pr, p = p, type = "Pr"))
}

p_list <- seq(30, 120, 10)

p1 <- ggplot(df_glm, aes(x = p-20, y = value, color = type)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = p_list) +
    scale_y_continuous(breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) +
    labs(x = "Noise Variables", y = "", color = "") +
    theme_classic() +
    theme(legend.position = "none")

p2 <- ggplot(df_glmnet, aes(x = p-20, y = value, color = type)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = p_list) +
    scale_y_continuous(breaks = seq(0, 1.0, 0.2),  limits = c(0, 1)) +
    labs(x = "Noise Variables", y = "", color = "") +
    theme_classic()

pdf("fig/example5.pdf", width = 5.7, height = 2)
p1 + p2
dev.off()