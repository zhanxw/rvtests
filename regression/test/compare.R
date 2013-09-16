files = commandArgs(trailingOnly = TRUE)[1:2]
options(stringsAsFactors = FALSE)

left <- read.table(files[1], header = FALSE)
right <- read.table(files[2], header = FALSE)

library(testthat)
expect_equal(left, right, tolerance = 1e-5)
