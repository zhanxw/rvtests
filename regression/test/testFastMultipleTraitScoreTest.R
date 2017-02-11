set.seed(12345)

N <- 10
M <- 1
C <- 2
T <- 3

cov <- matrix(rnorm(N*C), nr = N)
y <- matrix(rnorm(N*T), nr = N)
g <- matrix(rnorm(N*M), nr = N)

write.table(cov, "input.mt.cov", row.names = F, col.names = F)
write.table(y, "input.mt.y", row.names = F, col.names = F)
write.table(g, "input.mt.g", row.names = F, col.names = F)

##  full data set analysis
source("~/mylib/R/scoreTest-function.R")
idx <- complete.cases(cov[,1], y[,1], g)
linear.score(cov[idx,1], y[idx,1], g[idx,])

idx <- complete.cases(cov[,2], y[,2], g)
linear.score(cov[idx,2], y[idx,2], g[idx,])

idx <- complete.cases(cov, y[,2], g)
linear.score(cov[idx,], y[idx,2], g[idx,])

idx <- complete.cases(y[,1], g)
linear.score(NULL, y[idx,1], g[idx,])

## approx data set analysis
source("~/mylib/R/approxScoreTest-function.R")
debug(approx.linear.score)
approx.linear.score(cov[,1], y[,1], g)
approx.linear.score(cov[,2], y[,2], g)
approx.linear.score(cov, y[,2], g)
approx.linear.score(NULL, y[,1], g)

ret <- mapply(c, approx.linear.score(cov[,1], y[,1], g),
              approx.linear.score(cov[,2], y[,2], g),
              approx.linear.score(cov, y[,2], g),
              approx.linear.score(NULL, y[,1], g))
t(ret)
