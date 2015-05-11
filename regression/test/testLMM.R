param.true <- c(1, 1, 1)

genEst <- function(param) {
  set.seed(12345)

  n <- 200

  k1 <- matrix(runif(n*n, 0, 1), n, n)
  k1 <- cov2cor(t(k1) %*% k1)

  k2 <- matrix(runif(n*n, 0, 1), n, n)
  k2 <- cov2cor(t(k2) %*% k2)

  I <- diag(rep(1, n))

  Sigma <- param[1] * k1 + param[2] * k2 + param[3] * I
  summary(eigen(Sigma)$values)

  Cov <- cbind(rep(1,n), rnorm(100))
  beta <- c(1, 1)
  mu <- Cov %*% beta
  library(mvtnorm)
  y <- t(rmvnorm(1, mean = mu, sigma = Sigma))
  g <- matrix(sample(c(0, 1, 2), n, replace = TRUE), nc = 1)

  SigmaInv <- solve(Sigma)

  write.table(y, "input.lmm.y", quote = F, row.names = F, col.names = F)
  write.table(Cov, "input.lmm.x", quote = F, row.names = F, col.names = F)
  write.table(g, "input.lmm.g", quote = F, row.names = F, col.names = F)
  write.table(k1, "input.lmm.k1", quote = F, row.names = F, col.names = F)
  write.table(k2, "input.lmm.k2", quote = F, row.names = F, col.names = F)
}

getEst <- function(param) {
  x <-  as.matrix(read.table("input.lmm.x"))
  y <-  as.matrix(read.table("input.lmm.y"))
  g <-  as.matrix(read.table("input.lmm.g"))
  k1 <- as.matrix(read.table("input.lmm.k1"))
  k2 <- as.matrix(read.table("input.lmm.k2"))

  Cov <- x
  n <- nrow(x)
  I <- diag(rep(1, n))

  f <- function(param, mask) {
    print(param)
   print(mask)
    Sigma <- param[1] * k1 * mask[1] + param[2] * k2 * mask[2] + param[3] * I * mask[3]
    # Sigma <- param[1] * k1  + param[2] * k2  + param[3] * I * 0
    SigmaInv <- solve(Sigma)
    beta <- solve(t(Cov) %*% SigmaInv %*% Cov, t(Cov) %*% SigmaInv %*% y)
    resid <- y - Cov %*% beta
    ## llk <- -0.5 * (n*log(2*pi) + log(det(Sigma)) + t(resid) %*% SigmaInv %*% resid )
    ## llk <- -0.5 * (n*log(2*pi) + sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)) + t(resid) %*% SigmaInv %*% resid )
    L <- chol(Sigma)
    llk <- -0.5 * (n*log(2*pi) + sum(log(diag(L))) * 2 + t(resid) %*% SigmaInv %*% resid )
    llk <- -llk ## optim is a minimizer
    print(llk)
    return(llk)
  }

  #optim(param, f, method = "Nelder-Mead", lower = rep(1e-8, 3))
  optim(param, f, mask = c(1, 1, 0), method = "Nelder-Mead")
  optim(param, f, mask = c(1, 0, 1), method = "Nelder-Mead")

  getUV <- function(param) {
    print(param)
    Sigma <- param[1] * k1 + param[2] * k2 + param[3] * I
    SigmaInv <- solve(Sigma)
    beta <- solve(t(Cov) %*% SigmaInv %*% Cov, t(Cov) %*% SigmaInv %*% y)
    resid <- y - Cov %*% beta
    ## llk <- -0.5 * (n*log(2*pi) + log(det(Sigma)) + t(resid) %*% SigmaInv %*% resid )
    llk <- -0.5 * (n*log(2*pi) + sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)) + t(resid) %*% SigmaInv %*% resid )
    U <- t(g) %*% SigmaInv %*% (y - Cov %*% beta)
    V <- t(g) %*% SigmaInv %*% g - t(g) %*% SigmaInv %*% Cov %*% solve(t(Cov) %*% SigmaInv %*% Cov) %*% t(Cov) %*% SigmaInv %*% g
    print(sprintf("U = %g", U))
    print(sprintf("V = %g", V))
  }

  # 0.85413361 0.99349013 0.00000001
  # same
  optim(param, f, method = "L-BFGS-B", lower = rep(1e-8, 3))$par

  Sigma <- param[1] * k1 + param[2] * k2 + param[3] * I

  beta <- c(1, 1)
  mu <- Cov %*% beta
  library(mvtnorm)
  y <- t(rmvnorm(1, mean = mu, sigma = Sigma))
  g <- matrix(sample(c(0, 1, 2), n, replace = TRUE), nc = 1)

  SigmaInv <- solve(Sigma)

  U <- t(g) %*% SigmaInv %*% (y - Cov %*% beta)
  V <- t(g) %*% SigmaInv %*% g - t(g) %*% SigmaInv %*% Cov %*% solve(t(Cov) %*% SigmaInv %*% Cov) %*% t(Cov) %*% SigmaInv %*% g
  print(sprintf("U = %g", U))
  print(sprintf("V = %g", V))

}
# 0.85413361 0.99349013 0.00000001
# [1] 843.8605
getEst(param)

param.est <- getEst(param.true)

ret <- list()
for (i in seq(0.01, 100, length.out = 20)) {
  print(i)
  param.true <- c(1, 1, i)
  param.est <- getEst(param.true)
  ret[[length(ret) + 1]] <- c(param.true, param.est)
}

ret.all <- do.call(cbind, ret)
ret.all[c(1, 4), ]
ret.all[c(2, 5), ]
ret.all[c(3, 6), ]

ret.all.n100 <- ret.all
