require(Rcpp)
sourceCpp('methods/MAC.cpp')

MAC <- function(X, Y, X2, Y2, k = NULL) {
    if(is.vector(X)) {
        X <- as.matrix(X, ncol = 1)
        X2 <- as.matrix(X2, ncol = 1)
    }
    if(is.vector(Y)) {
        Y <- as.matrix(Y, ncol = 1)
        Y2 <- as.matrix(Y2, ncol = 1)
    }
    n <- dim(X)[1]
    m <- dim(X2)[2]
    if(is.null(k)) {
        k <- n
        list_i <- seq(1, n)
        list_j <- seq((n+1), (n + m))
    }else if(k <= n) {
        list_i <- order(X[1:n, 1])[round(seq(1, n, length.out = k))]
        list_j <- order(X2[1:m, 1])[round(seq(1, m, length.out = k))]
        list_j <- list_j + n
    }else {
        stop("k must be smaller than n")
    }
    mac <- compute_MAC(rbind(X, X2), rbind(Y, Y2), list_i, list_j, n)
    return(mac)
}

MAC2 <- function(x, y) {
    mac1 <- MAC(x[, 1], x[, 2], y[, 1], y[, 2])
    mac2 <- MAC(y[, 1], y[, 2], x[, 1], x[, 2])
    return(ifelse(mac1 > mac2, mac1, mac2))
}
