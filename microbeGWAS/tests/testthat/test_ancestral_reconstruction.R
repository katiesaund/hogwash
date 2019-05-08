library(microbeGWAS)
context("Ancestral reconstruction") -------------------------------------------#

# TODO I broke ancestral reconstruction by ML into several new subfunctions.
# Add unit tests for those subfunctions.

# test pick_recon_model --------------------------------------------------------
test_that("pick_recon_model gives an error when it claims to have a discrete phenotype, but continuous phenotype is given", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(1:num_cells, nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method))
})

test_that("pick_recon_model doesn't give an error when it is given a good discrete phenotype", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1,0), num_cells), nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method), NA)
})

test_that("pick_recon_model gives an error when it is given a continuous phenotype", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(1:num_cells, nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "continuous"
  index <- 1
  reconstruction_method <- "ML"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method))
})

test_that("pick_recon_model gives an error when it is a bad reconstruction method", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1,0), num_cells), nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "fake_method"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method))
})

test_that("pick_recon_model gives an error when it given a bad index", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1,0), num_cells), nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- num_col + 10
  reconstruction_method <- "ML"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method))
})

test_that("pick_recon_model gives an error when the input matrix is actually a vector", {
  num_col <- 9
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- rep(c(1, 0), num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  expect_error(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method))
})

test_that("pick_recon_model returns 'ER' for this test data", {
  num_col <- 9
  set.seed(1)
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1,0), num_cells), nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  expect_identical(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method), "ER")
})

test_that("pick_recon_model returns 'ARD' for this test particular data", {
  num_col <- 100
  set.seed(1)
  tree <- rtree(num_col, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1,1,1,1,1,1,1,0), num_cells), nrow = Ntip(tree), ncol = num_col)
  discrete_or_continuous <- "discrete"
  index <- 1
  reconstruction_method <- "ML"
  expect_identical(pick_recon_model(test_mat, tree, discrete_or_continuous, index, reconstruction_method), "ARD")
})

# test ancestral_reconstruction_by_ML ------------------------------------------
test_that("ancestral_reconstruction_by_ML with discrete input produce ancestral reconstruction with a value for each tip and node.", {
  set.seed(1)
  tree <- rtree(9, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_col <- 9
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rep(c(1, 0), num_cells), nrow = Ntip(tree), ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete")
  expected_length <- Ntip(tree) + Nnode(tree)
  expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
  expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
})

test_that("ancestral_reconstruction_by_ML with continuous input produce ancestral reconstruction with a value for each tip and node.", {
  set.seed(1)
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

test_that("ancestral_reconstruction_by_ML with continuous input produce ancestral reconstruction with a value for node.", {
  set.seed(1)
  tree <- rtree(9, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_col <- 9
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10), nrow = Ntip(tree), ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "continuous")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "continuous")
  expected_length <- Nnode(tree)
  expect_identical(length(dummy_pheno$node_anc_rec), expected_length)
  expect_identical(length(dummy_geno$node_anc_rec), expected_length)
})

test_that("ancestral_reconstruction_by_ML with continuous input produce ancestral reconstruction matrix with same dimensions as tree$edge matrix.", {
  set.seed(1)
  tree <- rtree(9, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_col <- 9
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10), nrow = Ntip(tree), ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "continuous")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "continuous")
  expected_rows <- nrow(tree$edge)
  expected_columns <- ncol(tree$edge)
  expect_identical(nrow(dummy_pheno$recon_edge_mat), expected_rows)
  expect_identical(ncol(dummy_geno$recon_edge_mat), expected_columns)
})

test_that("ancestral_reconstruction_by_ML with discrete input produces ancestral reconstruction with this known result.", {
  num_col <- 2
  num_tips <- 5
  set.seed(1)
  tree <- rtree(num_tips, rooted = TRUE)
  tree$node.label <- rep(100, Nnode(tree))
  num_cells <- num_col * Ntip(tree)
  test_mat <- matrix(c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1), nrow = Ntip(tree), ncol = num_col)
  dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete")
  dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete")
  expect_equivalent(dummy_pheno$tip_and_node_recon, c(1, 1, 1, 1, 0, 1, 1, 1, 1))
  expect_equivalent(dummy_pheno$node_anc_rec, c(1, 1, 1, 1))
  expect_equivalent(dummy_pheno$tip_and_node_rec_conf, c(1, 1, 1, 1, 1, 1, 1, 1, 0))
  expect_equivalent(dummy_pheno$recon_edge_mat[ , 1, drop = TRUE], c(1, 1, 1, 1, 1, 1, 1, 1))
  expect_equivalent(dummy_pheno$recon_edge_mat[ , 2, drop = TRUE], c(1, 1, 1, 1, 1, 1, 1, 0))
  expect_equivalent(dummy_geno$tip_and_node_recon,  c(0, 0, 0, 1, 1, 0, 0, 0, 1))
  expect_equivalent(dummy_geno$node_anc_rec, c(0, 0, 0, 1))
  expect_equivalent(dummy_geno$tip_and_node_rec_conf, c(1, 1, 1, 1, 1, 1, 1, 1, 0))
  expect_equivalent(dummy_geno$recon_edge_mat[ , 1, drop = TRUE], c(0, 0, 0, 0, 0, 0, 1, 1))
  expect_equivalent(dummy_geno$recon_edge_mat[ , 2, drop = TRUE], c(0, 0, 0, 0, 0, 1, 1, 1))
})

test_that("convert_to_edge_mat gives known result when given a valid tree and fake reconstruction", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  test_vector <- rep(1, Ntip(temp_tree) + Nnode(temp_tree))
  output_matrix <- matrix(1, ncol = 2, nrow = Nedge(temp_tree))
  expect_error(convert_to_edge_mat(temp_tree, test_vector), NA)
  expect_equivalent(convert_to_edge_mat(temp_tree, test_vector), output_matrix)


  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  test_vector <- c(1:sum(Ntip(temp_tree) + Nnode(temp_tree)))
  test_vector <- test_vector + 10
  output_matrix <- temp_tree$edge + 10
  expect_error(convert_to_edge_mat(temp_tree, test_vector), NA)
  expect_equivalent(convert_to_edge_mat(temp_tree, test_vector), output_matrix)
})

test_that("convert_to_edge_mat gives an error when not given a tree", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  test_vector <- rep(1, Ntip(temp_tree) + Nnode(temp_tree))
  output_matrix <- matrix(1, ncol = 2, nrow = Nedge(temp_tree))
  not_a_tree <- "foobar"
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- matrix(1, 10, 10)
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- NA
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
  not_a_tree <- NULL
  expect_error(convert_to_edge_mat(not_a_tree, test_vector))
})

test_that("convert_to_edge_mat gives an error when given vector or wrong dimension or type", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  not_the_reconstruction <- rep(1, Ntip(temp_tree))
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- NA
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- NULL
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- matrix(1, 10, 10)
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
  not_the_reconstruction <- "foobar"
  expect_error(convert_to_edge_mat(temp_tree, not_the_reconstruction))
})

# End script
