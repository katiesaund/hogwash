context("Genotype") #----------------------------------------------#

# test reduce_redundancy()

test_that("reduce_redundancy removes columns that have all ones or all zeroes, or all but one 1s.", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_mat <- matrix(c(1, 1, 1, 0, 0,
                       1, 1, 1, 1, 0,
                       1, 1, 1, 1, 1,
                       0, 0, 0, 0, 0,
                       0, 0, 1, 0, 0),
                     ncol = 5, nrow = 5)
  row.names(test_mat) <- temp_tree$tip.label
  colnames(test_mat) <- c("three_1_two_0", "four_1_one_0", "five_1_zero_0",
                          "zero_1_five_0", "one_1_four_0")
  temp_results <- reduce_redundancy(test_mat, temp_tree)
  expect_identical(temp_results$dropped_genotype_names, c("four_1_one_0", "five_1_zero_0", "zero_1_five_0", "one_1_four_0"))
  expect_identical(colnames(temp_results$mat), "three_1_two_0")
})

test_that("reduce_redundancy gives error is all columns contain only 0s.", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_mat <- matrix(0, ncol = 5, nrow = 5)
  row.names(test_mat) <- temp_tree$tip.label
  colnames(test_mat) <- letters[1:5]
  expect_error(reduce_redundancy(test_mat, temp_tree))
})

test_that("reduce_redundancy gives error is all columns contain letters.", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_mat <- matrix(letters[1:25], ncol = 5, nrow = 5)
  row.names(test_mat) <- temp_tree$tip.label
  colnames(test_mat) <- letters[1:5]
  expect_error(reduce_redundancy(test_mat, temp_tree))
})

test_that("reduce_redundancy removes no columns when all rows are two 1s and three 0s.", {
  set.seed(1)
  temp_tree <- ape::rtree(5)
  temp_tree$node.labels <- rep(100, ape::Nnode(temp_tree))
  test_mat <- matrix(c(0, 1, 1, 0, 0,
                       1, 1, 0, 0, 0,
                       1, 0, 0, 0, 1,
                       0, 0, 1, 1, 0,
                       0, 0, 0, 1, 1),
                     ncol = 5, nrow = 5)
  row.names(test_mat) <- temp_tree$tip.label
  colnames(test_mat) <- letters[1:5]
  temp_results <- reduce_redundancy(test_mat, temp_tree)
  expect_equal(length(temp_results$dropped_genotype_names), 0)
  expect_identical(colnames(temp_results$mat), letters[1:5])
})

