require(caret)
require(e1071)
require(randomForest)
require(ggplot2)
require(ggsci)
source("methods/simulation_methods.R")

# Example 3: crescent-shaped binary classification. The data are easy to
# classify with nonlinear models but structurally mismatched to a linear
# logistic regression boundary.
set.seed(123)
n <- 200
theta1 <- runif(n, -pi / 4, 3*pi/4)
x1 <- 5 * cos(theta1) + rnorm(n, sd = 0.5)
y1 <- 10 * sin(theta1) + rnorm(n, sd = 0.5)
class1 <- rep(0, n)

theta2 <- runif(n, -pi / 4, 3*pi/4)
x2 <- 5 * cos(theta2 + pi) + rnorm(n, 5, sd = 0.5)
y2 <- 10 * sin(theta2 + pi) + rnorm(n, -5, sd = 0.5)
class2 <- rep(1, n)

X <- rbind(cbind(x1, y1), cbind(x2, y2))
y <- c(class1, class2)
colnames(X) <- c("X1", "X2")

trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
trainData <- X[trainIndex, ]
trainLabels <- y[trainIndex]
testData <- X[-trainIndex, ]
testLabels <- y[-trainIndex]


# Working models with increasing structural flexibility.
lr_model <- glm(as.factor(y) ~ ., data = data.frame(X, y = as.factor(y)), family = "binomial")

# SVM with a radial kernel.
svm_rbf <- svm(as.factor(y) ~ ., data = data.frame(X, y = as.factor(y)), kernel = "radial", probability = TRUE)

# Random Forest baseline.
rf_model <- randomForest(as.factor(y) ~ ., data = data.frame(X, y = as.factor(y)))

# Predicted probabilities on the training subset for ps evaluation.
lr_train_prob <- predict(lr_model, type = "response", newdata = data.frame(trainData))

svm_rbf_train_prob <- attr(predict(svm_rbf, newdata = data.frame(trainData), probability = TRUE), "probabilities")[,2]

rf_train_prob <- predict(rf_model, type = "prob", newdata = data.frame(trainData))[,2]

# Test-set accuracy is reported together with ps to separate accuracy from
# structural optimality.
lr_pred <- predict(lr_model, newdata = data.frame(testData), type = "response")
lr_acc <- mean((lr_pred > 0.5) == testLabels)

svm_rbf_pred <- predict(svm_rbf, newdata = data.frame(testData))
svm_rbf_acc <- mean(svm_rbf_pred == testLabels)

rf_pred <- predict(rf_model, newdata = data.frame(testData))
rf_acc <- mean(rf_pred == testLabels)
cat("Accuracy on test set for linear model:", lr_acc, "\n")
cat("Accuracy on test set for SVM model:", svm_rbf_acc, "\n")
cat("Accuracy on test set for RF model:", rf_acc, "\n")
p_linar <- ps_with_probs(trainData, trainLabels, lr_train_prob)
p_svm <- ps_with_probs(trainData, trainLabels, svm_rbf_train_prob)
p_rf <- ps_with_probs(trainData, trainLabels, rf_train_prob)
cat("PS of Linear:", p_linar, "\n")
cat("PS of SVM:", p_svm, "\n")
cat("PS of RF:", p_rf, "\n")

# Build approximate decision-boundary points for plotting.
df <- data.frame(X, y = as.factor(y))
x1_grid <- seq(min(X[, 1]), max(X[, 1]), length.out = 100)
x2_grid <- seq(min(X[, 2]), max(X[, 2]), length.out = 100)
grid <- expand.grid(X1 = x1_grid, X2 = x2_grid)
grid$predicted <- predict(lr_model, type = "response", newdata = grid)
grid$predicted[order(abs(grid$predicted - 0.5))[1:150]]
grid <- grid[order(abs(grid$predicted - 0.5))[1:150], ]
grid$model <- rep("Linear Logistic regression", 150)
pred_df <- grid

grid <- expand.grid(X1 = x1_grid, X2 = x2_grid)
grid$predicted <- attr(predict(svm_rbf, newdata = data.frame(grid), probability = TRUE), "probabilities")[,2]
grid <- grid[order(abs(grid$predicted - 0.5))[1:150], ]
grid$model <- rep("Non-linear SVM", 150)
pred_df <- rbind(pred_df, grid)

grid <- expand.grid(X1 = x1_grid, X2 = x2_grid)
grid$predicted <- predict(rf_model, type = "prob", newdata = grid)[,2]
grid <- grid[order(abs(grid$predicted - 0.5))[1:150], ]
grid$model <- rep("Random Forest", 150)
pred_df <- rbind(pred_df, grid)

