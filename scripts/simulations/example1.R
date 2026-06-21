require(MASS)
require(glmnet)
require(ggplot2)
require(parallel)

# Example 1: show how Lasso feature recovery degrades as irrelevant
# variables are added to an otherwise linear signal.
n <- 200
p_base <- 20  # number of true signal variables
beta_true <- c(rep(1, p_base), rep(0, 500))
sigma <- 4
n_simulations <- 100  # number of repetitions

# Evaluate how many true variables are selected by cross-validated Lasso.
evaluate_lasso <- function(X, Y, beta_true) {
  # cv.glmnet selects the regularization parameter for Lasso.
  lasso_model <- cv.glmnet(X, Y, alpha = 1, standardize = TRUE)
  lasso_coef <- coef(lasso_model, s = "lambda.min")[-1]
  selected <- which(lasso_coef != 0)
  true_vars <- which(beta_true != 0)
  TP <- sum(selected %in% true_vars)
  return(TP)
}

# Number of noise variables added to the data-generating process.
noise_vars <- seq(50, 500, 50)

sub_test <- function(p1) {
  n <- 200
  p_base <- 20
  beta_true <- c(rep(1, p_base), rep(0, 500))
  sigma <- 4
  n_simulations <- 100
  p <- p_base + p1
  TP_simulations <- numeric(n_simulations)
  for (j in 1:n_simulations) {
    X <- matrix(rnorm(n * p), n, p)
    Y <- X[, 1:p_base] %*% beta_true[1:p_base] + rnorm(n, sd = sigma)
    TP_simulations[j] <- evaluate_lasso(X, Y, beta_true[1:p])
  }
  return(TP_simulations)
}

tps <- mclapply(noise_vars, FUN = sub_test, mc.cores = 10)
tps <- do.call(rbind, tps)

results_df <- data.frame(noise_vars = character(), TP = numeric())

for (i in seq_along(noise_vars)) {
  tmp_df <- data.frame(noise_vars = rep(noise_vars[i], 100), TP = tps[i, ])
  results_df <- rbind(results_df, tmp_df)
}

results_df$noise_vars <- factor(results_df$noise_vars, levels = noise_vars)
