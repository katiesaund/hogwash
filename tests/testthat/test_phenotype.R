context("Phenotype") #---------------------------------------------------------#
# test internal_report_phylogenetic_signal
test_that("internal_report_phylogenetic_signal gives no errors for correct inputs",  {
  # Discrete, random
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  row.names(pheno) <- tree$tip.label
  expect_error(internal_report_phylogenetic_signal(pheno, tree), NA)

  # discrete, BM
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2), nrow = num_tip, ncol = 1)
  pheno[pheno == "A", ] <- 0
  pheno[pheno == "B", ] <- 1
  storage.mode(pheno) <- "numeric"
  expect_error(internal_report_phylogenetic_signal(pheno, tree), NA)

  # continuous, random
  set.seed(1)
  pheno <- matrix(rnorm(n = num_tip, mean = 1, sd = 0.6))
  row.names(pheno) <- tree$tip.label
  expect_error(internal_report_phylogenetic_signal(pheno, tree), NA)

  # continuous, BM
  set.seed(1)
  pheno <- as.matrix(ape::rTraitCont(tree, model = "BM"), nrow = num_tip, ncol = 1)
  expect_error(internal_report_phylogenetic_signal(pheno, tree), NA)
})

test_that("internal_report_phylogenetic_signal gives errors for incorrect inputs",  {
  # no tree
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  set.seed(1)
  tree <- "foo"
  expect_error(internal_report_phylogenetic_signal(pheno, tree))

  # dataframe instead of matrix
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2),
                     nrow = num_tip,
                     ncol = 1)
  pheno[pheno == "A", ] <- 0
  pheno[pheno == "B", ] <- 1
  storage.mode(pheno) <- "numeric"
  pheno <- as.data.frame(pheno)
  expect_error(internal_report_phylogenetic_signal(pheno, tree))

  # characters instead of numerics
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2), nrow = num_tip, ncol = 1)
  expect_error(internal_report_phylogenetic_signal(pheno, tree))

  # vector instead of matrix
  set.seed(1)
  pheno <- ape::rTraitCont(tree, model = "BM")
  expect_error(internal_report_phylogenetic_signal(pheno, tree))
})

test_that("internal_report_phylogenetic_signal gives no errors for correct inputs, but not in the same order / unrooted",  {
  # not in the same order
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  row.names(pheno) <- paste0("t", 1:num_tip)
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  expect_error(internal_report_phylogenetic_signal(pheno, tree), NA)

})



## ------

# test report_phylogenetic_signal()
test_that("report_phylogenetic_signal gives no errors for correct inputs",  {
  # Discrete, random
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  row.names(pheno) <- tree$tip.label
  expect_error(report_phylogenetic_signal(pheno, tree), NA)

  # discrete, BM
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2), nrow = num_tip, ncol = 1)
  pheno[pheno == "A", ] <- 0
  pheno[pheno == "B", ] <- 1
  storage.mode(pheno) <- "numeric"
  expect_error(report_phylogenetic_signal(pheno, tree), NA)

  # continuous, random
  set.seed(1)
  pheno <- matrix(rnorm(n = num_tip, mean = 1, sd = 0.6))
  row.names(pheno) <- tree$tip.label
  expect_error(report_phylogenetic_signal(pheno, tree), NA)

  # continuous, BM
  set.seed(1)
  pheno <- as.matrix(ape::rTraitCont(tree, model = "BM"), nrow = num_tip, ncol = 1)
  expect_error(report_phylogenetic_signal(pheno, tree), NA)
})

test_that("report_phylogenetic_signal gives errors for incorrect inputs",  {
  # no tree
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  set.seed(1)
  tree <- "foo"
  expect_error(report_phylogenetic_signal(pheno, tree))

  # dataframe instead of matrix
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2),
                     nrow = num_tip,
                     ncol = 1)
  pheno[pheno == "A", ] <- 0
  pheno[pheno == "B", ] <- 1
  storage.mode(pheno) <- "numeric"
  pheno <- as.data.frame(pheno)
  expect_error(report_phylogenetic_signal(pheno, tree))

  # characters instead of numerics
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2), nrow = num_tip, ncol = 1)
  expect_error(report_phylogenetic_signal(pheno, tree))

  # vector instead of matrix
  set.seed(1)
  pheno <- ape::rTraitCont(tree, model = "BM")
  expect_error(report_phylogenetic_signal(pheno, tree))
})

test_that("report_phylogenetic_signal gives no errors for correct inputs, but not in the same order / unrooted",  {
  # not in the same order
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  row.names(pheno) <- paste0("t", 1:num_tip)
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  expect_error(report_phylogenetic_signal(pheno, tree), NA)

  # unrooted
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  row.names(pheno) <- paste0("t", 1:num_tip)
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  tree <- ape::unroot(tree)
  expect_error(report_phylogenetic_signal(pheno, tree), NA)
})


# calculate_lambda
test_that("calculate_lambda gives a lambda value for good data", {
  num_tip <- 25

  # continuous, random
  set.seed(1)
  pheno <- matrix(rnorm(n = num_tip, mean = 1, sd = 0.6))
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  row.names(pheno) <- tree$tip.label
  expect_error(calculate_lambda(pheno, tree), NA)
  expect_equal(round(calculate_lambda(pheno, tree), 5), .00007)

  # continuous, BM
  set.seed(1)
  pheno <- as.matrix(ape::rTraitCont(tree, model = "BM"), nrow = num_tip, ncol = 1)
  expect_error(calculate_lambda(pheno, tree), NA)
  expect_equal(round(calculate_lambda(pheno, tree), 5), 1.00064)
})

test_that("calculate_lambda gives error for bad data", {
  num_tip <- 25

  # phenotype / tree tips don't match
  set.seed(1)
  pheno <- matrix(rnorm(n = num_tip, mean = 1, sd = 0.6))
  row.names(pheno) <- paste0("t", 1:25)
  expect_error(calculate_lambda(pheno, tree))
})

# report_lambda
test_that("Each case of report_lambda is covered", {
  lambda <- -10
  expect_error(report_lambda(lambda), NA)

  lamdba <- 0
  expect_error(report_lambda(lambda), NA)

  lambda <- 0.5
  expect_error(report_lambda(lambda), NA)

  lambda <- 1
  expect_error(report_lambda(lambda), NA)

  lambda <- 10
  expect_error(report_lambda(lambda), NA)
})

test_that("Error message if not given a number", {
  lambda <- "foo"
  expect_error(report_lambda(lambda))

  lamdba <- matrix(0, 1, 1)
  expect_error(report_lambda(lambda))
})

# test calculate_d
test_that("calculate_d gives no errors when running on valid inputs", {
  # Discrete, random
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  row.names(pheno) <- tree$tip.label
  expect_error(calculate_d(pheno, tree), NA)
  set.seed(1)
  out <- round(calculate_d(pheno, tree), 2) # ~1.35
  expect_true(out < 1.5)
  expect_true(out > 1.2)

  # discrete, BM
  set.seed(210)
  pheno <- as.matrix(ape::rTraitDisc(tree, model = "ER", k = 2), nrow = num_tip, ncol = 1)
  pheno[pheno == "A", ] <- 0
  pheno[pheno == "B", ] <- 1
  storage.mode(pheno) <- "numeric"
  expect_error(calculate_d(pheno, tree), NA)
  set.seed(1)
  out <- round(calculate_d(pheno, tree), 2) # ~1.01
  expect_true(out < -.85)
  expect_true(out > -1.15)
})

test_that("calculate_d gives errors when running on invalid inputs", {
  # bad tree
  num_tip <- 25
  set.seed(1)
  pheno <- matrix(rbinom(num_tip, 1, 0.6))
  row.names(pheno) <- paste0("t", 1:num_tip)
  tree <- "foo"
  expect_error(calculate_d(pheno, tree))

  # bad pheno
  set.seed(1)
  tree <- ape::rcoal(n = num_tip)
  pheno <- "foo"
  expect_error(calculate_d(pheno, tree))
})


# report_d
test_that("Each case of report_d is covered", {
  dstat <- -10
  expect_error(report_d(dstat), NA)

  dstat <- 0
  expect_error(report_d(dstat), NA)

  dstat <- 0.5
  expect_error(report_d(dstat), NA)

  dstat <- 1
  expect_error(report_d(dstat), NA)

  dstat <- 10
  expect_error(report_d(dstat), NA)
})

test_that("Error message if not given a number", {
  dstat <- "foo"
  expect_error(report_d(dstat))

  dstat <- matrix(0, 1, 1)
  expect_error(report_d(dstat))
})

# test assign_pheno_type()
test_that("assign_pheno_type returns discrete for a discrete phenotype
          matrix", {
  phenotype <- matrix(0, nrow = 10, ncol = 1)
  expect_identical(assign_pheno_type(phenotype), "discrete")
})

test_that("assign_pheno_type returns continuous for a continuous phenotype
          matrix", {
  phenotype <- matrix(1.5, nrow = 10, ncol = 1)
  expect_identical(assign_pheno_type(phenotype), "continuous")
})

test_that("assign_pheno_type gives an error if phenotype is a vector, NA, or
          Null", {
  expect_error(assign_pheno_type(rep(1.5, 10)))
  expect_error(assign_pheno_type(NA))
  expect_error(assign_pheno_type(NULL))
})

test_that("assign_pheno_type gives an error if phenotype is a matrix of
          strings", {
  phenotype <- matrix(letters[1:10], nrow = 10, ncol = 1)
  expect_error(assign_pheno_type(phenotype))
})

test_that("assign_pheno_type gives an error if phenotype is a matrix of with two
          columns", {
  phenotype <- matrix(1:20, nrow = 10, ncol = 2)
  expect_error(assign_pheno_type(phenotype))
})

# test calculate_phenotype_change_on_edge
test_that("calculate_phenotype_change_on_edge returns a list of 5s, when first
          column is always greater than second column by 5", {
  temp_pheno <- matrix(1:10, ncol = 2)
  expect_equal(calculate_phenotype_change_on_edge(1:5, temp_pheno), rep(5, 5))
})

test_that("calculate_phenotype_change_on_edge returns a list of 5s, when first
          column is always less than second column by 5", {
  temp_pheno <- matrix(10:1, ncol = 2)
  expect_equal(calculate_phenotype_change_on_edge(1:5, temp_pheno), rep(5, 5))
})

test_that("calculate_phenotype_change_on_edge returns a list of 5s, when first
          column is always less than second column by 5", {
  temp_pheno <- matrix(-10:-1, ncol = 2)
  expect_equal(calculate_phenotype_change_on_edge(1:5, temp_pheno), rep(5, 5))
})

test_that("calculate_phenotype_change_on_edge returns a list of 5s, when first
          column is always less than second column by 5", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_equal(calculate_phenotype_change_on_edge(1:5, temp_pheno), rep(5, 5))
})

test_that("calculate_phenotype_change_on_edge returns a list of 5s, when first
column is always less than second column by 5, but
          this time only on a subset of edges", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_equal(calculate_phenotype_change_on_edge(2:4, temp_pheno), rep(5, 3))
})

test_that("calculate_phenotype_change_on_edge gives error for invalid input", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_error(calculate_phenotype_change_on_edge(9:10, temp_pheno))
})

test_that("calculate_phenotype_change_on_edge gives error for invalid input", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_error(calculate_phenotype_change_on_edge(matrix(NA), temp_pheno))
})

test_that("calculate_phenotype_change_on_edge gives error for invalid input", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_error(calculate_phenotype_change_on_edge(matrix(c(0, 0, 3)),
                                                  temp_pheno))
})



# test calc_raw_diff
test_that("calc_raw_diff returns a list of -5s, when first column is always
          greater than second column by 5", {
  temp_pheno <- matrix(1:10, ncol = 2)
  expect_equal(calc_raw_diff(1:5, temp_pheno), rep(-5, 5))
})

test_that("calc_raw_diff returns a list of 5s, when first column is always
          less than second column by 5", {
  temp_pheno <- matrix(10:1, ncol = 2)
  expect_equal(calc_raw_diff(1:5, temp_pheno), rep(5, 5))
})

test_that("calc_raw_diff returns a list of 5s, when first column is always less
          than second column by 5", {
  temp_pheno <- matrix(-10:-1, ncol = 2)
  expect_equal(calc_raw_diff(1:5, temp_pheno), rep(-5, 5))
})

test_that("calc_raw_diff returns a list of 5s, when first column is always less
          than second column by 5", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_equal(calc_raw_diff(1:5, temp_pheno), rep(-5, 5))
})

test_that("calc_raw_diff returns a list of 5s, when first column is always less
           than second column by 5, but this time only on a subset of edges", {
            temp_pheno <- matrix(-5:4, ncol = 2)
            expect_equal(calc_raw_diff(2:4, temp_pheno), rep(-5, 3))
})


test_that("calc_raw_diff gives error for invalid input", {
             temp_pheno <- matrix(-5:4, ncol = 2)
             expect_error(calc_raw_diff(10:12, temp_pheno))
})

test_that("calc_raw_diff gives error for invalid input", {
  temp_pheno <- matrix(-5:4, ncol = 2)
  expect_error(calc_raw_diff(matrix(c(1, 1, 1)), temp_pheno))
})

# # test check_if_convergence_occurs
# # check_if_convergence_occurs <- function(pheno, tree, continuous_or_discrete){
#
# test_that("check_if_convergence_occurs prints 'white noise model better' for
#            this specific test set", {
#   temp_pheno <- matrix(c(10, 1, 10, 1, 10, 1, 10, 1, 10, 1), ncol = 1)
#   set.seed(1)
#   temp_tree <- ape::rtree(10)
#   temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
#   temp_type <- "continuous"
#   row.names(temp_pheno) <- temp_tree$tip.label
#   expect_output(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type),
#                 "white noise model better")
# })
#
#
# test_that("check_if_convergence_occurs is silent for this specific test set", {
#   set.seed(1)
#   temp_tree <- ape::rtree(5)
#   temp_pheno <- matrix(rnorm(1:5), ncol = 1)
#   temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
#   temp_type <- "continuous"
#   row.names(temp_pheno) <- temp_tree$tip.label
#   expect_silent(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type))
# })
#
# test_that("check_if_convergence_occurs gives error when phenotype is
#            incorrectly formatted", {
#   set.seed(1)
#   temp_tree <- ape::rtree(5)
#   temp_pheno <- rnorm(1:5)
#   temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
#   temp_type <- "continuous"
#   expect_error(check_if_convergence_occurs(temp_pheno, temp_tree, temp_type))
# })
