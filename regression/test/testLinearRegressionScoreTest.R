outFile = commandArgs(trailingOnly = TRUE)[1]
setwd("~/rvtests/regression/test/")
source("ScoreTest.R")
set.seed(0)
n = 10
x = rnorm(n)
y = 1 +  0.5 * x + rnorm(n)

# output x, y 
X = cbind(rep(1,n), x)
write.table(file = "input.x", X, row.names = F, col.names =F)
write.table(file = "input.y", y, row.names = F, col.names =F)


conn = file(outFile, "w")

p.score = linear.score(Xcol=x, Y=y)$pvalue

for (i in 1:5){ # total 5 blocks in corresponding .cpp file
  cat("score_p\t",p.score,"\n", file = conn)
}

close(conn)

