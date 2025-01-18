require(ggplot2)
source("methods/simulation_methods.R")

#data generation
set.seed(100)
n <- 200
X <- runif(n, 0, 1)
e <- rnorm(n, 0, 1)
Y <- 10 * (X - 0.2)^2 + e

#true model
X_true <- seq(0, 1, length.out = 200)
Y_true <- 10 * (X_true - 0.2)^2

#linear model
linear_model <- lm(Y ~ X)


#model eval
summary(linear_model)
Residuals <- linear_model$residuals
Y2 <- linear_model$fitted.values
ps_linear <- ps_continuous(X, Y, Y2)
print(ps_linear)

quadratic_model <- lm(Y ~ X + I(X^2))
summary(quadratic_model)
Y3 <- quadratic_model$fitted.values
ps_quadratic <- ps_continuous(X, Y, Y3)
print(ps_quadratic)


#plot the result
p1 <- ggplot() +
  geom_point(aes(x = X, y = Y), size = 0.8) +
  geom_smooth(aes(x = X, y = Y), method = "lm", se = FALSE, colour = "blue") +
  geom_line(aes(x = X_true, y = Y_true), colour = "red", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

pdf("fig/example1.pdf", height = 2.5, width = 3)
p1
dev.off()
