setwd("~/rvtests/regression/test/")
source("ScoreTest.R")
set.seed(0)
n = 30
x = rnorm(n)
cov = rnorm(n)
y = 1 +  0.2 * x + cov + rnorm(n)

# output x, y 
X = cbind(rep(1,n), x)
write.table(file = "input.x", X, row.names = F, col.names =F)
write.table(file = "input.cov", cov, row.names = F, col.names =F)
write.table(file = "input.y", y, row.names = F, col.names =F)

# wald
ret = lm(y~x+cov)
summary(ret)
beta = coef(ret)
v = vcov(ret)
p.wald = coef(summary(ret))[2,4]

conn = file("output.R.lm", "w")

cat("wald_beta\t", file = conn)
cat(beta, file = conn, append = TRUE)
cat("\n", file = conn)

cat("wald_vcov\t", file = conn)
cat(v, file = conn, append = TRUE)
cat("\n", file = conn)

cat("wald_p\t", file = conn)
cat(p.wald, file = conn, append = TRUE)
cat("\n", file = conn)

#permutation
beta.real = coef(lm(y~x+cov))[2]
x.resid = residuals( lm(x~cov) )
y.resid = residuals( lm(y~cov) )

permutated_beta<-function(a){
  y.sample = sample(y.resid)
  coef(lm(y.sample~x.resid))[2]
}

beta.perm = sapply(seq(1000), permutated_beta)

p.perm = sum(abs(beta.perm) >= abs(beta.real)) / length(beta.perm) 
p.perm
cat("permutation_p\t", file = conn)
cat(p.perm, file = conn, append = TRUE)
cat("\n", file = conn)

p.score = linear.score(Cov=cov, Xcol=x, Y=y)$pvalue
p.score
cat("score_p\t", file = conn)
cat(p.score, file = conn, append = TRUE)
cat("\n", file = conn)

close(conn)
