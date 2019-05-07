library(microbeGWAS)
context("Phenotype") #----------------------------------------------#


# test assign_pheno_type()
test_that("assign_pheno_type returns discrete for a discrete phenotype matrix", {
  phenotype <- matrix(0, nrow = 10, ncol = 1)
  expect_identical(assign_pheno_type(phenotype), "discrete")
})

test_that("assign_pheno_type returns continuous for a continuous phenotype matrix", {
  phenotype <- matrix(1.5, nrow = 10, ncol = 1)
  expect_identical(assign_pheno_type(phenotype), "continuous")
})

test_that("assign_pheno_type gives an error if phenotype is a vector, NA, or Null", {
  expect_error(assign_pheno_type(rep(1.5, 10)))
  expect_error(assign_pheno_type(NA))
  expect_error(assign_pheno_type(NULL))
})

test_that("assign_pheno_type gives an error if phenotype is a matrix of strings", {
  phenotype <- matrix(letters[1:10], nrow = 10, ncol = 1)
  expect_error(assign_pheno_type(phenotype))
})

test_that("assign_pheno_type gives an error if phenotype is a matrix of with two columns", {
  phenotype <- matrix(1:20, nrow = 10, ncol = 2)
  expect_error(assign_pheno_type(phenotype))
})

# test calculate_phenotype_change_on_edge
# calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){

# test calc_raw_diff
# calc_raw_diff <- function(edge_list, ph_edges){

# test check_if_phenotype_normal
test_that("check_if_phenotype_normal prints a statement that the data is not normal when the input is 1:100", {
  phenotype <- matrix(1:100, ncol = 1)
  expect_output(check_if_phenotype_normal(phenotype, "continuous"))
})

test_that("check_if_phenotype_normal prints a statement when the input is all the same number", {
  phenotype <- matrix(10, ncol = 1)
  expect_error(check_if_phenotype_normal(phenotype, "continuous"))
})

test_that("check_if_phenotype_normal doesn't give any message when phenotype is from a normal distribution", {
  set.seed(1)
  phenotype <- matrix(rnorm(n = 100), ncol = 1)
  expect_silent(check_if_phenotype_normal(phenotype, "continuous"))
})

test_that("check_if_phenotype_normal does nothing when phenotype is 'discrete'", {
  phenotype <- matrix(c(1, 0, 1, 1, 1, 1, 1, 0), ncol = 1)
  expect_silent(check_if_phenotype_normal(phenotype, "discrete"))
})

# test check_if_convergence_occurs
# check_if_convergence_occurs <- function(pheno, tree, continuous_or_discrete){

test_that("check_if_convergence_occurs prints 'white noise model better' for this specific test set", {
  temp_pheno <- matrix(c(10, 1, 10, 1, 10, 1, 10, 1, 10, 1), ncol = 1)
  set.seed(1)
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_type <- "continuous"
  row.names(temp_pheno) <- temp_tree$tip.label
  expect_output(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type), "white noise model better")
})


test_that("check_if_convergence_occurs is silent for this specific test set", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_pheno <- matrix(rnorm(1:5), ncol = 1)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_type <- "continuous"
  row.names(temp_pheno) <- temp_tree$tip.label
  expect_silent(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type))
})

test_that("check_if_convergence_occurs gives error when phenotype is incorrectly formatted", {
  set.seed(1)
  temp_tree <- rtree(5)
  temp_pheno <- rnorm(1:5)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_type <- "continuous"
  expect_error(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type))
})

#TODO add checks for when inputs are bad to check_if_convergence_occurs()
