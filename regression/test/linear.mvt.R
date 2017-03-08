#' VT test for continous trait
#' sample size: n
#' 
#' @param y.res, vector of length n, residual = (y - y.null)
#' @param S, n by p matrix, genotype matrix
#' @param covar, n by d matrix, covariate matrix (need to include intercept)
#' @param v, vector of length 1, variance is the same in linear model
#' @return a list
linear.mvt <- function(y.res, S, covar, v) {
  if (length(v) > 1) {
    stopifnot(all(v == v[1]))
  }
    
  if (FALSE) {
    U = c(sum(y.res * S[,1]), sum(y.res * S[,2]) )
    U
  }
  U <- y.res %*% S  ## faster
  # print(U)
  U
  
  if (FALSE) {  ## naiive calculation of vsz
    tmp <- NULL
    for( i in 1:N) {
      if (i == 1) {
        tmp <- covar[i, ] * v[i] * S[i, 1]
      } else {
        tmp <- tmp + covar[i, ] * v[i] * S[i, 1]
      }
    }
  }
  #vsz.1 <- t(covar) %*% (v*S[,1]) ## faster
  #vsz.2 <- t(covar) %*% (v*S[,2]) ## faster
  ##vsz <- t(covar) %*% diag(v) %*% S ## fastest = cbind(vsz.1, vsz.2)
  vsz <- t(covar) %*% S ## *v ## fastest = cbind(vsz.1, vsz.2)
  
  zzt <- t(covar) %*% covar ## * v 
  if (FALSE) {
    V <- c( (sum(v * S[,1]^2) - t(vsz.1) %*% solve(zzt) %*% (vsz.1)),
            (sum(v * S[,2]^2) - t(vsz.2) %*% solve(zzt) %*% (vsz.2)))  
  }
  V <- v * (colSums((S*S)) - colSums(vsz * (solve(zzt) %*% vsz)))  # fastest
  V
  
  if (FALSE) {
    U.1 <- y.res * (S[,1] - covar %*% vsz.1 / zzt)
    U.2 <- y.res * (S[,2] - covar %*% vsz.2 / zzt)
    V.all <- matrix( c(sum(U.1 * U.1), sum(U.1 * U.2), sum(U.1*U.2), sum(U.2*U.2)), 2)
    V.all ## faster
  }
  
  Uk  <- diag(y.res) %*% (S - covar %*% solve(zzt) %*% vsz) ## Uk <- cbind(U.1, U.2)
  V.all <- t(Uk) %*% Uk ## fastest
  
  t.all <- U /sqrt(V)
  t.max <- max(abs(U/sqrt(V)))
  t.max
  
  library(mvtnorm)
  #print(t.max)
  #print(S)
  p <- pmvnorm(lower = -rep(t.max, ncol(S)), 
               upper = rep(t.max, ncol(S)), 
               mean = rep(0, ncol(S)),
               sigma = (V.all)
  )
  return (list(U = U, V = V, t = t.all, t.max = t.max, V.all = V.all, p = 1 - p))
}
