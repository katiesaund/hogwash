library(microbeGWAS)
context("Discrete algorithm") -------------------------------------------------#
# TODO finishing writing unit tests
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
  num_genotypes <- 15
  set.seed(1)
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_trans <- temp_hi_conf_edges <- rep(list(0), num_genotypes)
  for (k in 1:num_genotypes){
    temp_geno_trans[[k]] <- c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
    temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  temp_geno_trans[[15]] <- c(0, 0, 1, 1, 0, 0, 0, 1, 0, 0)
  temp_pheno_trans <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  temp_geno <- matrix(1, ncol = num_genotypes, nrow = Ntip(temp_tree)) # doesn't match recon or transition, just made up for now.
  temp_perm <- 8
  temp_fdr <- 0.25

  expect_error(discrete_calculate_pvals(temp_geno_trans, temp_pheno_trans, temp_tree, temp_geno, temp_perm, temp_fdr, temp_hi_conf_edges), NA)

  disc_trans_results <- discrete_calculate_pvals(temp_geno_trans, temp_pheno_trans, temp_tree, temp_geno, temp_perm, temp_fdr, temp_hi_conf_edges)

  expect_equal(round(as.numeric(disc_trans_results$hit_pvals[1]), 3), 0.444)
  expect_equal(disc_trans_results$observed_overlap[1], 1)
  expect_equal(length(disc_trans_results$permuted_count[[1]]), temp_perm)
})


# test discrete_permutation

test_that("discrete_permutation returns expected results given this dummy data", {
  num_samples <- 6
  num_genotypes <- 15
  set.seed(1)
  temp_tree <- rtree(num_samples)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  num_edge <- Nedge(temp_tree)
  temp_geno_trans <- temp_hi_conf_edges <- rep(list(0), num_genotypes)
  for (k in 1:num_genotypes){
    temp_geno_trans[[k]] <- c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
    temp_hi_conf_edges[[k]] <- rep(1, num_edge)
  }
  temp_geno_trans[[15]] <- c(0, 0, 1, 1, 0, 0, 0, 1, 0, 0)
  temp_pheno_trans <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  temp_geno <- matrix(1, ncol = num_genotypes, nrow = Ntip(temp_tree)) # doesn't match recon or transition, just made up for now.
  temp_perm <- 8

  temp_num_edges_with_geno_trans <- sapply(temp_geno_trans, function(x) sum(x))
  temp_num_hi_conf_edges <- sapply(temp_hi_conf_edges, function(x) sum(x))

  for (m in 1:num_genotypes){
    permuted_geno <- discrete_permutation(temp_tree, temp_perm, temp_num_edges_with_geno_trans, temp_num_hi_conf_edges, num_edge, temp_hi_conf_edges, m)
    expect_equal(nrow(permuted_geno), temp_perm)
    expect_equal(ncol(permuted_geno), num_edge)
    expect_lt(max(rowSums(permuted_geno)), temp_num_edges_with_geno_trans[m] + 1)
  }

  edges_selected <- NULL
  for (r in 1:nrow(permuted_geno)){
    temp <-   c(1:10)[permuted_geno[r, ] == 1]
    edges_selected <- c(edges_selected, temp)
  }

  edge_probability <- temp_tree$edge.length / sum(temp_tree$edge.length) # high confidence is one for all genotypes, so can do just one calculation for all genotypes in this case.
  edges_selected_other_way <- sample(1:num_edge, size = length(edges_selected), replace = TRUE, prob = edge_probability)
  ks_results <-  withCallingHandlers(ks.test(edges_selected_other_way, edges_selected),
                                     warning=function(w) {
                                       if (grepl("cannot compute exact p-value with ties", w$message))
                                         invokeRestart("muffleWarning")
                                     } )
  expect_true(ks_results$p.value > 0.05)

})

# test count_empirical_both_present

test_that("count_empirical_both_present gives X given Y", {
  temp_pheno_vec <- c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
  temp_hi_conf_edge <- NULL
  temp_hi_conf_edge[[1]] <- c(0, 1 ,1, 0, 0, 1, 1, 0, 0, 1, 1, 1)
  temp_hi_conf_edge[[2]] <- c(0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1)
  temp_permuted_mat <- matrix(rbinom(120, 1, .5), ncol = 12, nrow = 10)
  temp_index <- 1
  expect_error(count_empirical_both_present(temp_permuted_mat, temp_pheno_vec, temp_hi_conf_edge, temp_index), NA)
  expect_warning(count_empirical_both_present(temp_permuted_mat, temp_pheno_vec, temp_hi_conf_edge, temp_index), NA)
})

# test count_empirical_only_geno_present

# count_empirical_only_geno_present <- function(permuted_mat, emp_both_present){
#   # Function description -------------------------------------------------------
#   # Find number of edges on the tree where the permuted genotype is a transition
#   # but the phenotype is not (either a transition phenotype edge or phenotype
#   # reconstruction present -- it depends on the type of test being run).
#   #
#   # Inputs:
#   # permuted_mat. Matrix. Nrow = number of permutations. Ncol = number of edges. Either 0 or 1.
#   # emp_both_present. Numeric vector.
#   #
#   # Outputs:
#   # result. Numeric vector. Length = num perm.
#   #
#   # Check input ----------------------------------------------------------------
#   check_if_binary_matrix(permuted_mat)
#   if (length(emp_both_present) != nrow(permuted_mat)){
#     stop("Input dimension mismatch")
#   }
#
#   # Function -------------------------------------------------------------------
#   result <- sapply(1:nrow(permuted_mat), function(x) {
#     sum(permuted_mat[x, ]) - emp_both_present[x]
#   })
#   # Check & return output --------------------------------------------------------------
#   if (nrow(permuted_mat) != length(result)){
#     stop("result dimension mismatch")
#   }
#
#   return(result)
# } # end count_empirical_only_geno_present()
#
