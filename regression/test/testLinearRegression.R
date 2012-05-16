setwd("~/rvtests/regression/test/")
source("ScoreTest.R")
set.seed(0)
n = 30
x = rnorm(n)
y = 1 +  0.5 * x + rnorm(n)

# output x, y 
X = cbind(rep(1,n), x)
write.table(file = "input.x", X, row.names = F, col.names =F)
write.table(file = "input.y", y, row.names = F, col.names =F)

# wald
ret = lm(y~x)
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
beta.real = coef(lm(y~x))[2]

permutated_beta<-function(a){
  y.sample = sample(y)
  coef(lm(y.sample~x))[2]
}

beta.perm = sapply(seq(200), permutated_beta)

p.perm = sum(abs(beta.perm) >= abs(beta.real)) / length(beta.perm) 
p.perm
cat("permutation_p\t", file = conn)
cat(p.perm, file = conn, append = TRUE)
cat("\n", file = conn)

p.score = linear.score(Xcol=x, Y=y)$pvalue
cat("score_p\t", file = conn)
cat(p.score, file = conn, append = TRUE)
cat("\n", file = conn)

close(conn)
