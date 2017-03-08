outFile = commandArgs(trailingOnly = TRUE)[1]

source("logistic.mvt.R")
set.seed(0)
N <- 1000
M <- 2
p <- seq(M) * 0.01

x <- do.call(cbind, lapply(p, function(x){rbinom(N, 2, x)}))
colMeans(x)
covar <- matrix(rnorm(N * 3), N, 3)
covar <- cbind(rep(1, nrow(covar)), covar)
y <- covar %*% seq(ncol(covar))
y <- 1/(1+exp(-y))
y <- rbinom(N, 1, y)

summary(glm(y ~ x + covar - 1, family = "binomial"))
colSums(x)

y.null <- predict(glm(y ~ covar - 1, family = "binomial"), type = "response")
summary(y.null)
v <- y.null * (1-y.null)
y.res <- y - y.null

S <- x
S[,2] <- S[,1] + S[,2]
x <- S

ret <- logistic.mvt(y.res, S, covar, v)

out <- rbind(ret$U, ret$V, ret$t, ret$V.all, c(ret$p, 0))

write.table(x,     "input.logistic.mvt.x", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(y,     "input.logistic.mvt.y", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(covar, "input.logistic.mvt.cov", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(out, outFile, quote = FALSE, row.names = FALSE, col.names=FALSE)

if (F){
    l.res <-as.matrix(read.table("res"))
    l.S <-  as.matrix(read.table("S"))
}
