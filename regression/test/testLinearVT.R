outFile = commandArgs(trailingOnly = TRUE)[1]
library(mvtnorm)
source("linear.mvt.R")
set.seed(0)
N <- 1000
M <- 2
p <- seq(M) * 0.01

x <- do.call(cbind, lapply(p, function(x){rbinom(N, 2, x)}))
colMeans(x)
cor.rand <- matrix(rnorm(3*3), 3)
cor.rand <- t(cor.rand) %*% cor.rand
covar <- rmvnorm(N, mean = rep(0, 3), sigma = cor.rand)
covar <- cbind(rep(1, nrow(covar)), covar)
#y <- covar %*% seq(ncol(covar)) + rnorm(N)
y <- rnorm(N)

summary(ret <- lm(y ~ covar - 1))
colSums(x)

y.res <- resid(ret)
v <- sum(y.res ^ 2) / N

S <- x
S[,2] <- S[,1] + S[,2]
x <- S

ret <- linear.mvt(y.res, S, covar, v)

out <- rbind(ret$U, ret$V, ret$t, ret$V.all, c(ret$p, 0))

write.table(x,     "input.linear.mvt.x", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(y,     "input.linear.mvt.y", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(covar, "input.linear.mvt.cov", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(out, outFile, quote = FALSE, row.names = FALSE, col.names=FALSE)

if (F){
    l.res <-as.matrix(read.table("res"))
    l.S <-  as.matrix(read.table("S"))
}
