# Helper functions for computing the structure-improvement p-value (ps)
# under continuous outcomes, GLMs, glmnet fits, or supplied probabilities.
source("methods/NewMac.R")
require(parallel)

# Continuous response: bootstrap residuals from the working model and
# compare observed data against model-generated data with MAC.
BootStrapTb <- function(Predictor, TrueValue, Estimator, BootNumber) {
    eta <- TrueValue - Estimator
    if(is.vector(eta)) {
        n <- length(eta)
    }else {
        n <- dim(eta)[1]
    }
    if(n > 2000) {
        k <- 200
    }else {
        k <- max(c(100, n))
    }
    boot_tb <- function(boot, Predictor, TrueValue, Estimator, eta, n, k) {
        set.seed(boot + 100)
        eta_star <- sample(eta, n, replace = TRUE)
        Y2 <- Estimator + eta_star
        set.seed(boot + 99)
        eta_star2 <- sample(eta, n, replace = TRUE)
        Y22 <- Estimator + eta_star2
        mac1 <- MAC(Predictor, TrueValue, Predictor, Y2, k)
        mac2 <- MAC(Predictor, Y22, Predictor, Y2, k)
        return(c(mac1, mac2))
    }
    Tb <- mclapply(1:BootNumber, boot_tb, Predictor = Predictor,
        TrueValue = TrueValue, Estimator = Estimator, eta = eta, n = n, k = k, mc.cores = 10)
    Tb <- do.call(rbind, Tb)
    return(Tb)
}
ps_continuous <- function(X, Y, Y2, BootNumber = 1000) {
  Tb <- BootStrapTb(X, Y, Y2, BootNumber)
  T0 <- mean(Tb[, 1])
  null_t <- Tb[, 2]
  P <- (sum(null_t >= T0) + 1) / (BootNumber + 1)
  return(P)
}


# Binary GLM: parametric bootstrap from fitted probabilities.
ps_glm <- function(X, Y, model, BootNumber = 1000) {
  predicted_probs <- predict(model, type = "response")
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


# glmnet logistic model: extract the selected lambda and bootstrap labels.
compute_P_glmnet <- function(X, Y, fit, BootNumber = 1000) {
  beta_hat <- fit$glmnet.fit$beta[, fit$index[1]]
  predicted_probs <- 1 / (1 + exp(- X %*% beta_hat - fit$glmnet.fit$a0[fit$index[1]]))
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

# Generic binary model interface when predicted probabilities are already
# available from an external model.
ps_with_probs <- function(X, Y, predicted_probs, BootNumber = 1000) {
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
