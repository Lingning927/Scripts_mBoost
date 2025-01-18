library(data.table)

X_all <- fread("scripts/simulations/example3/X.csv", header = TRUE)
Y_all <- fread("scripts/simulations/example3/Y.csv", header = TRUE)

set.seed(100)
test_id <- sample(seq(1, nrow(X_all)), round(0.2*nrow(X_all)))

X_train <- X_all[-test_id, ]
Y_train <- Y_all[-test_id, ]
write.csv(X_train, "scripts/simulations/example3/X_train.csv", row.names = FALSE)
write.csv(Y_train, "scripts/simulations/example3/Y_train.csv", row.names = FALSE)

X_test <- X_all[test_id, ]
Y_test <- Y_all[test_id, ]
write.csv(X_test, "scripts/simulations/example3/X_test.csv", row.names = FALSE)
write.csv(Y_test, "scripts/simulations/example3/Y_test.csv", row.names = FALSE)

X_train <- fread("scripts/simulations/example3/X_train.csv", header = TRUE)
