options(stringsAsFactors = FALSE)
default.args <- c("testBoltLMM.out", "correct.testBoltLMM.out")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- default.args
}
# print(args)
stopifnot(all(file.exists(args)))

f1 <- read.table(args[1], header = TRUE)
f2 <- read.table(args[2], header = TRUE)
stopifnot(dim(f1)  == dim(f2))

min.corr <- min(sapply(seq(ncol(f1)), function(x) { cor(f1[,x], f2[,x])}))
print(sprintf("min correlation = %g", min.corr))
stopifnot(min.corr > 0.95)



