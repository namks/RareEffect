## Data generation------------------------------------------------------------------------------------------------------------------------
set.seed(0)

n <- 2000   # sample size
p <- 100    # number of SNP in a gene


## ---------------------------------------------------------------------------------------------------------------------------------------
X <- matrix(rbinom(n * 2, 1, 0.5), nrow = n) # Matrix of covariates
alpha <- c(1, 1)   # effect size of covariates


## ---------------------------------------------------------------------------------------------------------------------------------------
G <- matrix(0, nrow = n, ncol = p)  # Matrix of genotypes
AF <- runif(p, min = 0, max = 1)    # AF of each SNP

for (j in 1:p) {
    G[,j] = rbinom(n, size = 2, prob = AF[j])
}


## ---------------------------------------------------------------------------------------------------------------------------------------
tau <- 1    # variance component of random term
K <- diag(1, nrow = p)  # beta ~ MVN(0, tau * K)
beta <- t(chol(tau * K)) %*% rnorm(p)   # random effect terms


## ---------------------------------------------------------------------------------------------------------------------------------------
Y <- X %*% alpha + G %*% beta + rnorm(n, 0, 1)


## ---------------------------------------------------------------------------------------------------------------------------------------
#library(Matrix)
GKG <- G %*% K %*% t(G)


## EMMAX from SKAT------------------------------------------------------------------------------------------------------------------------
library(SKAT)
system.time({out<-SKAT_NULL_emmaX (Y ~ X, K=GKG)})


## ---------------------------------------------------------------------------------------------------------------------------------------
print(out$ve / out$va)
print(out$va)
print(out$ve)


## EMMAX ---------------------------------------------------------------------------------------------------------------------------------
ll <- function(delta) {
    H <- GKG + diag(delta, n)
    d <- determinant(H, logarithm = T)
    
    Hinv_X <- solve(H, X)
    Hinv_Y <- solve(H, Y)
    Hinv_e <- (Hinv_Y - Hinv_X %*% solve(t(X) %*% Hinv_X) %*% t(X) %*% Hinv_Y)
    
    R <- t(Y - X %*% solve(t(X) %*% Hinv_X) %*% t(X) %*% Hinv_Y) %*% Hinv_e
    sigma_sq <- R / n
    
    LL <- -0.5 * (n * log(2 * pi * sigma_sq) + d$modulus + n)
    return (as.numeric(LL))
}

sigma_sq <- function(delta){
    H <- GKG + diag(delta, n)
    
    Hinv_X <- solve(H, X)
    Hinv_Y <- solve(H, Y)
    Hinv_e <- (Hinv_Y - Hinv_X %*% solve(t(X) %*% Hinv_X) %*% t(X) %*% Hinv_Y)
    R <- t(Y - X %*% solve(t(X) %*% Hinv_X) %*% t(X) %*% Hinv_Y) %*% Hinv_e
    sigma_sq <- R / n
    return (as.numeric(sigma_sq))
}

opt <- optimize(ll, lower=0, upper=10, maximum=TRUE)
opt_delta <- opt$maximum
va <- sigma_sq(opt_delta)
ve <- opt_delta * va


## ---------------------------------------------------------------------------------------------------------------------------------------
print(opt_delta)
print(va)
print(ve)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Using optimal delta, calculate alpha = (X' H^{-1} X)^{-1} X' H^{-1} y
alpha.hat <- solve(t(X) %*% solve(GKG + diag(opt_delta, n)) %*% X) %*% t(X) %*% solve(GKG + diag(opt_delta, n)) %*% Y
print(alpha.hat)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Using optimal delta, calculate the posterior beta given y
beta.post <- K %*% t(G) %*% solve(GKG + diag(opt_delta, n)) %*% (Y - X %*% alpha.hat)
print(beta.post)


## FaST-LMM ------------------------------------------------------------------------------------------------------------------------------
# SVD of W
L <- chol(K)    # S = LL^T
W <- G %*% L    # GKG^T = G L L^T G^T = W W^T 
svd <- svd(W)
U <- svd$u
S <- (svd$d)^2
V <- svd$v

UtY <- t(U) %*% Y
UtX <- t(U) %*% X


## ---------------------------------------------------------------------------------------------------------------------------------------
# Calculate likelihood
ll <- function(delta) {
    a_hat1 <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    for (i in 1:p) {
        a_hat1 <- a_hat1 + (UtX[i,] %*% t(UtX[i,])) / (S[i] + delta)
    }
    for (i in 1:n) {
        a_hat1 <- a_hat1 + (1 / delta) * (X - U %*% UtX)[i,] %*% t((X - U %*% UtX)[i,])
    }
    
    a_hat2 <- matrix(0, nrow = ncol(X), ncol = 1)
    for (i in 1:p) {
        a_hat2 <- a_hat2 + (UtX[i,] %*% t(UtY[i,])) / (S[i] + delta)
    }
    for (i in 1:n) {
        a_hat2 <- a_hat2 + (1 / delta) * (X - U %*% UtX)[i,] %*% t((Y - U %*% UtY)[i,])
    }
    
    a_hat <- solve(a_hat1) %*% a_hat2
    
    Lik1 <- 0
    for (i in 1:p) {
        Lik1 <- Lik1 + (UtY - UtX %*% a_hat)[i,]^2 / (S[i] + delta)
    }
    Lik1 <- Lik1 + 1 / delta * sum((Y - U %*% UtY - (X - U %*% UtX) %*% a_hat)^2)

    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - p) * log(delta) + n + n * log(1 / n * Lik1))
    
    return (as.numeric(out))
}

a.hat <- function(delta) {
    a_hat1 <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    for (i in 1:p) {
        a_hat1 <- a_hat1 + (UtX[i,] %*% t(UtX[i,])) / (S[i] + delta)
    }
    for (i in 1:n) {
        a_hat1 <- a_hat1 + (1 / delta) * (X - U %*% UtX)[i,] %*% t((X - U %*% UtX)[i,])
    }
    
    a_hat2 <- matrix(0, nrow = ncol(X), ncol = 1)
    for (i in 1:p) {
        a_hat2 <- a_hat2 + (UtX[i,] %*% t(UtY[i,])) / (S[i] + delta)
    }
    for (i in 1:n) {
        a_hat2 <- a_hat2 + (1 / delta) * (X - U %*% UtX)[i,] %*% t((Y - U %*% UtY)[i,])
    }
    
    a_hat <- solve(a_hat1) %*% a_hat2
    
    return (as.numeric(a_hat))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Estimate beta

## Find optimal delta (using likelihood)
opt <- optimize(ll, lower=0, upper=10, maximum=TRUE)
opt_delta <- opt$maximum

## beta = [sum 1/S+ii + delta (U_1^T X) ~~]


## ---------------------------------------------------------------------------------------------------------------------------------------
# Using optimal delta, calculate alpha
alpha.hat <- a.hat(opt_delta)
print(alpha.hat)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Using optimal delta, calculate the posterior beta given y
beta.post <- K %*% t(G) %*% solve(GKG + diag(opt_delta, n)) %*% (Y - X %*% alpha.hat)
print(beta.post)


## Moment-based approach -----------------------------------------------------------------------------------------------------------------
resid1 <- resid(lm(Y ~ X))

len.e <- length(resid1) * (length(resid1) - 1) / 2
e.ij <- rep(0, len.e)
k.ij <- rep(0, len.e)
idx <- 1

for (i in 1:(length(resid1)-1)) {
    for (j in (i+1):length(resid1)) {
        e.ij[idx] <- resid1[i] * resid1[j]
        k.ij[idx] <- GKG[i, j]
        idx <- idx + 1
    }
}


## ---------------------------------------------------------------------------------------------------------------------------------------
m1 <- lm(e.ij ~ k.ij)
summary(m1)


## ---------------------------------------------------------------------------------------------------------------------------------------
m2 <- lm(e.ij ~ k.ij - 1)
summary(m2)


