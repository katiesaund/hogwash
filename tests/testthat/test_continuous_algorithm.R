context("Continuous algorithm") #----------------------------------------------#

# test get_sig_hits_while_correcting_for_multiple_testing
test_that("get_sig_hit_and_mult_test_corr works with valid input", {
  set.seed(1)
  fake_p <- rnorm(n = 10, mean = 0.2, sd = 0.2)
  fake_fdr <- 0.1
  results <- get_sig_hit_and_mult_test_corr(fake_p, fake_fdr)
  expect_identical(round(results$hit_pvals$fdr_corrected_pvals, 7),
                   c(0.2490308,
                     0.3862944,
                     0.1795316,
                     0.5190562,
                     0.3862944,
                     0.1795316,
                     0.3862944,
                     0.3862944,
                     0.3862944,
                     0.3473058))
  expect_equal(nrow(results$sig_pvals), 0)
})

# test calc_sig
test_that("calc_sig gives expected results given valid inputs", {
  num_isolates <- 5
  num_loci <- 8
  set.seed(1)
  temp_tree <- ape::rtree(num_isolates)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_geno <- matrix(c(0, 1), nrow = num_isolates, ncol = num_loci)
  temp_perm <- 100
  temp_geno_trans <- temp_conf <- temp_geno_recon <- rep(list(NULL), num_loci)
  for (i in 1:num_loci) {
    temp_geno_trans[[i]]$transition <- c(1, 0, 1, 0, 1, 0, 1, 0)
    temp_geno_trans[[i]]$trans_dir <- c(1, 0, -1, 0, 1, 0, -1, 0)
    temp_conf[[i]] <- c(1, 1, 1, 0, 0, 1, 1, 1)
  }
  set.seed(1)
  temp_pheno <- matrix(rnorm(ape::Nedge(temp_tree) * 2),
                       ncol = 2,
                       nrow = ape::Nedge(temp_tree))

  hi_conf_list <- NULL
  hi_conf_list$genotype <- temp_geno
  hi_conf_list$high_conf_ordered_by_edges <- temp_conf
  hi_conf_list$genotype_transition <- temp_geno_trans
  hi_conf_list$tr_and_pheno_hi_conf <- rep(1, ape::Nedge(temp_tree))

  results <- calc_sig(hi_conf_list,
                      temp_perm,
                      temp_tree,
                      temp_pheno)
  expect_equal(results$num_genotypes, num_loci)
  expect_equal(length(results$observed_pheno_non_trans_delta[[1]]), 4)
  expect_equal(length(results$observed_pheno_trans_delta[[1]]), 4)
})


test_that("calc_sig returns a lower p-value when phenotype change is
           noticeably higher on transition edges than on non-transition tree
          edges", {
  num_isolates <- 5
  num_loci <- 2
  set.seed(1)
  temp_tree <- ape::rtree(num_isolates)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_geno <- matrix(c(0, 1), nrow = num_isolates, ncol = num_loci)
  temp_perm <- 100
  temp_geno_trans <- temp_conf <- rep(list(NULL), num_loci)
  for (i in 1:num_loci) {
    temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1, 1, 0)
    temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1, 1, 1)
    temp_conf[[i]]                  <- c(1, 1, 1, 1, 1, 1, 1, 1)
  }
  set.seed(1)
  temp_pheno <-
    matrix(c(-1, -1, 10, -.5, 10, 10, 0, 0, 0, -.5, 11, 30, 100, 40, 150, 0.1),
           ncol = 2,
           nrow = ape::Nedge(temp_tree))
  # Transition edges deltas: 30, 30.5, 90, 150
  # Non-transition edge deltas: 0.5, 1, 1, 0.1

  hi_conf_list <- NULL
  hi_conf_list$genotype <- temp_geno
  hi_conf_list$high_conf_ordered_by_edges <- temp_conf
  hi_conf_list$genotype_transition <- temp_geno_trans
  hi_conf_list$tr_and_pheno_hi_conf <- rep(1, ape::Nedge(temp_tree))

  results <- calc_sig(hi_conf_list,
                      temp_perm,
                      temp_tree,
                      temp_pheno)
  alpha <- 0.15
  expect_true(results$pvals[1] < alpha)
  expect_equal(results$observed_pheno_non_trans_delta[[1]],
               abs(temp_pheno[c(1:3, 8), 1] - temp_pheno[c(1:3, 8), 2]))
  expect_equal(round(results$observed_pheno_trans_delta[[1]], 3),
               abs(temp_pheno[4:7, 1] - temp_pheno[4:7, 2]))
})

test_that("calc_sig returns a non-significant p-value when phenotype change is
          identical on transition and non-transition edges", {
  num_isolates <- 4
  num_loci <- 2
  set.seed(1)
  temp_tree <- ape::rtree(num_isolates)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_geno <- matrix(c(0, 1), nrow = num_isolates, ncol = num_loci)
  temp_perm <- 100
  temp_geno_trans <- temp_conf <- rep(list(NULL), num_loci)
  for (i in 1:num_loci) {
    temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1)
    temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1)
    temp_conf[[i]]                  <- c(1, 1, 1, 1, 1, 1)
  }
  set.seed(1)
  temp_pheno <- matrix(c(-15, -15, 10, 0, 10, 10, -10, 0, 7, 5, 25, 7),
                       ncol = 2,
                       nrow = ape::Nedge(temp_tree))

  hi_conf_list <- NULL
  hi_conf_list$genotype <- temp_geno
  hi_conf_list$high_conf_ordered_by_edges <- temp_conf
  hi_conf_list$genotype_transition <- temp_geno_trans
  hi_conf_list$tr_and_pheno_hi_conf <- rep(1, ape::Nedge(temp_tree))

  # Transition edges deltas: 5 15 3
  # Non-transition edge deltas: 5 15 3
  results <- calc_sig(hi_conf_list,
                      temp_perm,
                      temp_tree,
                      temp_pheno)
  alpha <- 0.01
  expect_true(results$pvals[1] > alpha)

  expect_equal(results$observed_pheno_non_trans_delta[[1]],
               abs(temp_pheno[1:3, 1] - temp_pheno[1:3, 2]))
  expect_equal(round(results$observed_pheno_trans_delta[[1]], 3),
               abs(temp_pheno[4:6, 1] - temp_pheno[4:6, 2]))
})

test_that("calc_sig returns an error when the only confident edges are
          transition edges. (All non-transition edges are low confidence)", {
  num_isolates <- 4
  num_loci <- 2
  set.seed(1)
  temp_tree <- ape::rtree(num_isolates)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))
  temp_geno <- matrix(c(0, 1), nrow = num_isolates, ncol = num_loci)
  temp_perm <- 100
  temp_geno_trans <- temp_conf <- rep(list(NULL), num_loci)
  for (i in 1:num_loci) {
    temp_geno_trans[[i]]$transition <- c(0, 0, 0, 1, 1, 1)
    temp_geno_trans[[i]]$trans_dir  <- c(1, 1, 1, 1, 1, 1)
    temp_conf[[i]]                  <- c(0, 0, 0, 1, 1, 1)
  }
  set.seed(1)
  temp_pheno <- matrix(c(-15, -15, 10, 0, 10, 10, -10, 0, 7, 5, 25, 7),
                       ncol = 2,
                       nrow = ape::Nedge(temp_tree))

  hi_conf_list <- NULL
  hi_conf_list$genotype <- temp_geno
  hi_conf_list$high_conf_ordered_by_edges <- temp_conf
  hi_conf_list$genotype_transition <- temp_geno_trans
  hi_conf_list$tr_and_pheno_hi_conf <- rep(1, ape::Nedge(temp_tree))
  # Transition edges deltas: 5 15 3
  # Non-transition edge deltas: 5 15 3
  expect_error(calc_sig(hi_conf_list,
                        temp_perm,
                        temp_tree,
                        temp_pheno))
})


# test continuous_permutation
test_that("continuous_permutation gives valid results", {
  num_isolates <- 40
  num_loci <- 8
  set.seed(1)
  temp_tree <- ape::rtree(num_isolates)
  temp_tree$node.label <- rep(100, ape::Nnode(temp_tree))

  # Make each edge hi conf (tree, pheno, and geno)
  # Assign the delta to each edge as the edge number (for simplicity)
  temp_conf_index <- pheno_delta_list <- non_trans_index <- rep(list(NULL), num_loci)
  for (i in 1:num_loci) {
    temp_conf_index[[i]] <- pheno_delta_list[[i]] <- c(1:ape::Nedge(temp_tree))
    non_trans_index[[i]] <- c(1:40)
  }

  num <- 1
  perm <- 100
  results <- continuous_permutation(non_trans_index, temp_tree, temp_conf_index, perm, num, pheno_delta_list)
  odd_numbers <- seq(from = 1, to = 77, by = 2)
  expect_equal(length(results), perm)
  expect_true(max(results) <= ape::Nedge(temp_tree))
  expect_true(min(results) >= 1)
})
