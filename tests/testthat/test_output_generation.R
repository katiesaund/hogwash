context("Output generation") #-------------------------------------------------#

test_that("create_contingency_table() doesn't give error for synchronous" ,{
  geno_edges <- list(c(0, 0, 0, 0), c(0, 0, 1, 1), c(1, 1, 1, 1))
  pheno_edges <- c(1, 0, 1, 0)
  geno <- matrix(0, ncol = 3)
  colnames(geno) <- letters[1:3]

  sync_names <- list(c("geno_transition", "geno_not_trans"),
                     c("pheno_transition", "pheno_not_trans"))
  expected_results <- rep(list(NULL), length(geno_edges))
  expected_results[[1]] <-
    matrix(c(0, 2, 0, 2), nrow = 2, ncol = 2, dimnames = sync_names)
  expected_results[[2]] <-
    matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2, dimnames = sync_names)
  expected_results[[3]] <-
    matrix(c(2, 0, 2, 0), nrow = 2, ncol = 2, dimnames = sync_names)
  names(expected_results) <- letters[1:3]

  temp_results <- create_contingency_table(geno_edges,
                                           pheno_edges,
                                          geno,
                                          "synchronous")

  expect_equal(temp_results, expected_results)
})

test_that("create_contingency_table() doesn't give error for phyc" ,{
  geno_edges <- list(c(0, 0, 0, 0), c(0, 0, 1, 1), c(1, 1, 1, 1))
  pheno_edges <- c(1, 0, 1, 0)
  geno <- matrix(0, ncol = 3)
  colnames(geno) <- letters[1:3]
  phyc_names <- list(c("geno_transition", "geno_not_trans"),
                     c("pheno_present", "pheno_absent"))
  expected_results <- rep(list(NULL), length(geno_edges))
  expected_results[[1]] <-
    matrix(c(0, 2, 0, 2), nrow = 2, ncol = 2, dimnames = phyc_names)
  expected_results[[2]] <-
    matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2, dimnames = phyc_names)
  expected_results[[3]] <-
    matrix(c(2, 0, 2, 0), nrow = 2, ncol = 2, dimnames = phyc_names)
  names(expected_results) <- letters[1:3]

  temp_results <- create_contingency_table(geno_edges,
                                           pheno_edges,
                                           geno,
                                           "phyc")
  expect_equal(temp_results, expected_results)
})
