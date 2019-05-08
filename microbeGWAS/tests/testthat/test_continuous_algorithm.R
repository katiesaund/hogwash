library(microbeGWAS)
context("Continuous algorithm") -----------------------------------------------#

test_that("run_ks_test returns a known test statistic for given data.", {
  temp_pheno <- temp_pheno <- matrix(rnorm(20), ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:8)
  non_trans_index <- c(9:10)
  result <- run_ks_test(trans_index, non_trans_index, temp_pheno)
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno), NA)
  expect_equivalent(0.625, result$statistic)
})

test_that("run_ks_test gives an error when the indices are too few to run ks.test()", {
  temp_pheno <- temp_pheno <- matrix(1:20, ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:9)
  non_trans_index <- ""
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno))
})
