# check different regressions
logistic.score <- function(Cov=NULL, Y, Xcol){
  if (is.null(Cov)){
    Xnull = matrix(rep(1, length(Y)), ncol = 1)
  } else {
    Xnull = cbind(rep(1,nrow(as.matrix(Cov))), Cov)
  }
  Z = t(Xnull)
  S = t(Xcol)  
  n = nrow(Xnull)
  d = ncol(Xnull)
  b.null = coef(glm(Y~Xnull-1, family=binomial))
  p.null = 1/(1+exp(- Xnull %*% b.null))
  v.null = (exp(- Xnull %*% b.null))/ (1+exp(- Xnull %*% b.null))^2
  U = sum( (Y - p.null) * Xcol)
  V = sum ( v.null*Xcol*Xcol )
  left = matrix(0, nrow = d, ncol = 1)
  mid = matrix(0, nrow = d, ncol = d)
  for (i in 1:nrow(Xnull)){
    Zi = matrix(Xnull[i,], ncol = 1)
    left = left + v.null[i] * S[i] * Zi
    mid = mid + v.null[i] * Zi %*% t(Zi)
  }
  V = V - t(left) %*% solve(mid) %*% left
  pvalue = 1 - pchisq( U*U/V, df=1)
  list(U=U, V=V, pvalue = pvalue)
}

linear.score <- function(Cov=NULL, Y, Xcol){
  if (is.null(Cov)){
    Xnull = matrix(rep(1, length(Y)), ncol = 1)
  } else {
    Xnull = cbind(rep(1,nrow(as.matrix(Cov))), Cov)
  }
  Z = t(Xnull)
  S = t(Xcol)  
  n = nrow(Xnull)
  d = ncol(Xnull)
  b.null = coef(lm(Y~Xnull-1))
  p.null = Xnull %*% b.null
  v.null = sum( (Y - p.null)^2) / n
  U = sum( (Y - p.null) * Xcol)
  V = sum ( Xcol*Xcol )
  left = matrix(0, nrow = d, ncol = 1)
  mid = matrix(0, nrow = d, ncol = d)
  for (i in 1:nrow(Xnull)){
    Zi = matrix(Xnull[i,], ncol = 1)
    left = left + S[i] * Zi
    mid = mid + Zi %*% t(Zi)
  }
  V = v.null * (V - t(left) %*% solve(mid) %*% left)
  pvalue = 1 - pchisq( U*U/V, df=1)
  list(U=U, V=V, pvalue = pvalue, sigma2 = v.null)  
}

# logistic regression 
set.seed(0)
y = rep(c(0,1), c(5,5))
x = rep(0,10)
x[2] = -1
x[3] = 0.5
x[7] = 2
x[8] = 6
x[9] = 1.5

z = rep(c(1,2,3,4,5), 2)

# wald
summary(glm(y~x, family = binomial))

# score
logistic.score(Y = y, Xcol = x)

# wald + cov
summary(glm(y~x+z, family = binomial))

# score + cov
logistic.score(Cov = z, Y = y, Xcol = x)

# linear regression 
# wald
summary(lm(y~x))

# score
linear.score(Y = y, Xcol = x)

# wald + cov
summary(lm(y~x+z))

# score + cov
linear.score(Cov = z, Y = y, Xcol = x)


ret = glm(y~z)