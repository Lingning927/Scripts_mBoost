require(MASS)
require(glmnet)
require(ggplot2)
require(parallel)

n <- 200
p_base <- 20  #true related variables
beta_true <- c(rep(1, p_base), rep(0, 500))
sigma <- 4
n_simulations <- 100  #number of repetitions

#evaluate the feature selection ability of lasso
evaluate_lasso <- function(X, Y, beta_true) {
  #cv.glmnet for Lasso
  lasso_model <- cv.glmnet(X, Y, alpha = 1, standardize = TRUE)
  lasso_coef <- coef(lasso_model, s = "lambda.min")[-1]
  selected <- which(lasso_coef != 0)
  true_vars <- which(beta_true != 0)
  TP <- sum(selected %in% true_vars)
  return(TP)
}

#the number of nosie variables
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

pdf("fig/example4.pdf", height = 3, width = 5)
ggplot(results_df, aes(x = noise_vars, y = TP)) +
  geom_boxplot(fill = "#0975c2b9") +
  scale_y_continuous(limits = c(10, 20)) +
  labs(x = "Number of Noisy Variables",
       y = "True Positives") +
  theme_classic()
dev.off()