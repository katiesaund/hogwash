library(microbeGWAS)

# TODO
# This barely covers any test cases for this function. Need to expand tests.

context("ancestral_reconstruction_by_ML") #-------------------------------------------------#
test_that("ancestral_reconstruction_by_ML with discrete input produce ancestral reconstruction with a value for each tip and node.", {
  tree <- rtree(9, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_col <- 9
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells), nrow = Ntip(tree), ncol = num_col)
  check_if_binary_matrix(test_mat)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete")
  expected_length <- Ntip(tree) + Nnode(tree)
  expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
})

test_that("ancestral_reconstruction_by_ML with continuous input produce ancestral reconstruction with a value for each tip and node.", {
  tree <- rtree(9, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_col <- 9
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10), nrow = Ntip(tree), ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "continuous")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "continuous")
  expected_length <- Ntip(tree) + Nnode(tree)
  expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
})
