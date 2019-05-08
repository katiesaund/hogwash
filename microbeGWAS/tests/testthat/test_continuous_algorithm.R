library(microbeGWAS)
context("Continuous algorithm") -----------------------------------------------#

# test run_ks_test
test_that("run_ks_test returns a known test statistic for given data.", {
  set.seed(1)
  temp_pheno <- temp_pheno <- matrix(rnorm(20), ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:8)
  non_trans_index <- c(9:10)
  answer <- 0.375
  names(answer) <- "D"
  result <- run_ks_test(trans_index, non_trans_index, temp_pheno)
  expect_equivalent(answer, result$statistic)
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno), NA)
})

test_that("run_ks_test gives an error when the indices are too few to run ks.test()", {
  temp_pheno <- temp_pheno <- matrix(1:20, ncol = 2)
  num_edges <- nrow(temp_pheno)
  trans_index <- c(1:9)
  non_trans_index <- ""
  expect_error(run_ks_test(trans_index, non_trans_index, temp_pheno))
})

# test get_sig_hits_while_correcting_for_multiple_testing
test_that("get_sig_hits_while_correcting_for_multiple_testing gives known adjusted p-values for given data", {
  set.seed(1)
  fake_p <- rnorm(n = 10, mean = 0.2, sd = 0.2)
  fake_fdr <- 0.1
  results <- get_sig_hits_while_correcting_for_multiple_testing(fake_p, fake_fdr)
  expect_identical(round(results$hit_pvals$fdr_corrected_pvals, 7), c(0.2490308, 0.3862944, 0.1795316, 0.5190562, 0.3862944, 0.1795316, 0.3862944, 0.3862944, 0.3862944, 0.3473058))
  expect_equal(nrow(results$sig_pvals), 0)
})
