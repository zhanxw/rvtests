outFile = commandArgs(trailingOnly = TRUE)[1]

library(logistf)
data(sex2)

y <- sex2$case
x <- cbind(1, subset(sex2, select = c("age", "oc", "vic", "vicl", "vis", "dia")))
write.table(x, "input.firth.x", quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(y, "input.firth.y", quote = FALSE, row.names = FALSE, col.names=FALSE)

fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2, pl = FALSE)
#str(fit)
#coef(fit)

## x2 <- as.matrix(read.table("input.firth.x", header = FALSE))
## x2 <- x2[,-1]
## y2 <- as.matrix(read.table("input.firth.y", header = FALSE))
## coef(logistf(y2 ~ x2 ))

##

beta<-coef(fit)
beta.se<-sqrt(diag(vcov(fit)))
#str(summary(fit))
out <- cbind(beta, beta.se)

write.table(out, outFile, quote = FALSE, row.names = FALSE, col.names=FALSE)

## ## test
if (FALSE) {
    d <- as.matrix(read.table("matD"))
    covB <- as.matrix(read.table("matCovB"))
    diag(d %*% covB)
    norm( d %*% covB )
}
