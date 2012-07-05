outFile = commandArgs(trailingOnly = TRUE)[1]
setwd("~/rvtests/regression/test/")
source("ScoreTest.R")
set.seed(0)
n = 200
x = rnorm(n)
y = rbinom(n, 1, 1/ ( 1 + exp( - (1 +  0.5 *x))))

# output x, y 
X = cbind(rep(1,n), x)
write.table(file = "input.x", X, row.names = F, col.names =F)
write.table(file = "input.y", y, row.names = F, col.names =F)


conn = file(outFile, "w")

p.score = logistic.score(Xcol=x, Y=y)$pvalue
for (i in 1:5){ # total 5 blocks in corresponding .cpp file
  cat("score_p\t",p.score,"\n", file = conn)
}

close(conn)

