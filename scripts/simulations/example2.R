require(ggplot2)
source("methods/simulation_methods.R")

# Example 2: polynomial regression where a linear working model is
# intentionally misspecified and a quadratic model matches the truth.
set.seed(100)
n <- 200
X <- runif(n, 0, 1)
e <- rnorm(n, 0, 1)
Y <- 10 * (X - 0.2)^2 + e

# Ground-truth curve used for plotting.
X_true <- seq(0, 1, length.out = 200)
Y_true <- 10 * (X_true - 0.2)^2

# Fit and evaluate the linear working model.
linear_model <- lm(Y ~ X)


# ps <= 0.05 suggests that the model structure remains improvable.
summary(linear_model)
Residuals <- linear_model$residuals
Y2 <- linear_model$fitted.values
ps_linear <- ps_continuous(X, Y, Y2)
print(ps_linear)

# Fit and evaluate the correctly specified quadratic model.
quadratic_model <- lm(Y ~ X + I(X^2))
summary(quadratic_model)
Y3 <- quadratic_model$fitted.values
ps_quadratic <- ps_continuous(X, Y, Y3)
print(ps_quadratic)


# Plot observations, the linear fit, and the true curve.
p1 <- ggplot() +
  geom_point(aes(x = X, y = Y), size = 0.8) +
  geom_smooth(aes(x = X, y = Y), method = "lm", se = FALSE, colour = "blue") +
  geom_line(aes(x = X_true, y = Y_true), colour = "red", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))
