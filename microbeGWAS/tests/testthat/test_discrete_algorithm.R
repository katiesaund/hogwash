library(microbeGWAS)
context("Discrete algorithm") -------------------------------------------------#

# test calculate_permutation_based_p_value
test_that("calculate_permutation_based_p_value returns a significant p-value when statistics are much lower than observed", {
  nperm <- 1000
  permuted_tests <- rep(0.1, nperm)
  real_test <- 0.9
  alpha <- 0.01
  expect_true(calculate_permutation_based_p_value(permuted_tests, real_test, nperm) < alpha)
})

test_that("calculate_permutation_based_p_value returns a non-significant p-value when statistics are cenetered around observed", {
  nperm <- 1000
  permuted_tests <- rnorm(n = nperm, mean = 0)
  real_test <- 0
  alpha <- 0.01
  expect_true(calculate_permutation_based_p_value(permuted_tests, real_test, nperm) > alpha)
})

# test count_hits_on_edges
test_that("count_hits_on_edges returns 3 edges shared and 7 edges only with genotype given this test data", {
  num_samples <- 6
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_recon <- temp_hi_conf_edges <- rep(list(0), num_samples)
  for (k in 1:num_samples){
    temp_geno_recon[[k]] <- temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  num_pheno_and_geno_present <- 3
  num_just_geno_present <- num_edge - num_pheno_and_geno_present
  temp_pheno_recon <- c(rep(1, num_pheno_and_geno_present), rep(0, num_just_geno_present))

  results <- count_hits_on_edges(temp_geno_recon, temp_pheno_recon, temp_hi_conf_edges, temp_tree)
  expect_equal(results$both_present[1], num_pheno_and_geno_present)
  expect_equal(results$only_geno_present[1], num_just_geno_present)
})

test_that("count_hits_on_edges returns 0 edges shared and 0 edges only with genotype given this all absent genotype", {
  num_samples <- 6
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_recon <- temp_hi_conf_edges <- rep(list(0), num_samples)
  for (k in 1:num_samples){
    temp_geno_recon[[k]] <- temp_hi_conf_edges[[k]] <- rep(0, num_edge)
  }
  temp_pheno_recon <- rep(1, num_edge)
  results <- count_hits_on_edges(temp_geno_recon, temp_pheno_recon, temp_hi_conf_edges, temp_tree)
  expect_equal(results$both_present[1], 0)
  expect_equal(results$only_geno_present[1], 0)
})


# test discrete_calculate_pvals
test_that("discrete_calculate_pvals returns expected results given this dummy data", {
  num_samples <- 6
  num_genotypes <- 10
  set.seed(1)
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_trans <- temp_hi_conf_edges <- rep(list(0), num_genotypes)
  for (k in 1:num_genotypes){
    temp_geno_trans[[k]] <- c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
    temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  temp_geno_trans[[10]] <- c(0, 0, 1, 1, 0, 0, 0, 1, 0, 0)
  temp_pheno_trans <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  temp_geno <- matrix(1, ncol = num_genotypes, nrow = Ntip(temp_tree)) # doesn't match recon or transition, just made up for now.
  temp_perm <- 8
  temp_fdr <- 0.25
  disc_trans_results <- discrete_calculate_pvals(temp_geno_trans, temp_pheno_trans, temp_tree, temp_geno, temp_perm, temp_fdr, temp_hi_conf_edges)

  expect_equal(round(as.numeric(disc_trans_results$hit_pvals[1]), 3), 0.763)
  expect_equal(disc_trans_results$observed_overlap[1], 1)
  expect_equal(length(disc_trans_results$permuted_count[[1]]), temp_perm)
})
