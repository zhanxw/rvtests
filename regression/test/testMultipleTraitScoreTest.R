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

source("~/mylib/R/scoreTest-function.R")


scale <- function(g){
  apply(g, 2, function(x) {
  a <- x - mean(x)
  a / sqrt(sum(a^2))
})}
cov <- scale(cov)
y <- scale(y)
g <- scale(g)

linear.score(cov[,1], y[,1], g)
linear.score(cov[,2], y[,2], g)
linear.score(cov, y[,2], g)
linear.score(NULL, y[,1], g)

sum(y[,1]*g)
sum(cov[,1]*g)
sum(cov[,1]*y[,1])
