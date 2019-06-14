library(microbeGWAS)
context("Validation") #--------------------------------------------------------#

# test check_dimensions --------------------------------------------------------
test_that("check_dimensions gives an error for when given a dataframe instead of a matrix", {
  # Set up
  temp_mat <- as.data.frame(matrix(1:100, nrow = 10, ncol = 10))

  # Test
  expect_error(check_dimensions(temp_mat, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given a vector instead of a matrix", {
  # Set up
  temp_vector <- c(1:100)

  # Test
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given min_rows is too large", {
  # Set up
  num_row <- 10
  temp_mat <- matrix(rep("A", 100), nrow = num_row, ncol = 10)

  # Test
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = 10,
                                min_rows = num_row + 99,
                                exact_cols = NULL,
                                min_cols = 1))
})

test_that("check_dimensions gives an error for when given a min_cols it too large", {
  # Set up
  num_col <- 10
  temp_mat <- matrix(rep("A", 100), nrow = 10, ncol = num_col)

  # Test
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = NULL,
                                min_rows = 1,
                                exact_cols = NULL,
                                min_cols = num_col + 2))
})

test_that("check_dimensions gives an error for when given exact_cols that is smaller or larger than ncol()", {
  # Set up
  num_col <- 10
  temp_mat <- matrix(rep("A", 100), nrow = 2, ncol = num_col)

  # Test
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = NULL,
                                min_rows = 1,
                                exact_cols = num_col - 1,
                                min_cols = 1))
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = NULL,
                                min_rows = 1,
                                exact_cols = num_col + 1,
                                min_cols = 1))
})

test_that("check_dimensions gives an error for when given exact_rows that is smaller or larger than nrow()", {
  # Set up
  num_row <- 10
  temp_mat <- matrix(rep("A", 100), nrow = num_row, ncol = 1)

  # Test
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = num_row - 1,
                                min_rows = 1,
                                min_cols = 1))
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = num_row + 1,
                                min_rows = 1,
                                min_cols = 1))
})

test_that("check_dimensions gives an error for when given NULL instead of a matrix", {
  # Set up
  temp_vector <- NULL

  # Test
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given NA instead of a matrix", {
  # Set up
  temp_vector <- NA

  # Test
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions doesn't give an error for when given a valid matrix", {
  # Set up
  temp_matrix <- matrix(0, nrow = 10, ncol = 5)

  # Test
  expect_error(check_dimensions(temp_vector, 10, 1, 5, 1))
})

# test check_num_between_0_and_1 -----------------------------------------------
test_that("check_num_between_0_and_1 gives errors for negative numbers", {
  # Test
  expect_error(check_num_between_0_and_1(-1))
  expect_error(check_num_between_0_and_1(-Inf))
  expect_error(check_num_between_0_and_1(-.0000001))
})

test_that("check_num_between_0_and_1 does not give errors for valid inputs", {
  # Test
  expect_error(check_num_between_0_and_1(0), NA)
  expect_error(check_num_between_0_and_1(1), NA)
  expect_error(check_num_between_0_and_1(1/2), NA)
  expect_error(check_num_between_0_and_1(0.00000001), NA)
})

test_that("check_num_between_0_and_1 gives errors for numbers larger than 1", {
  # Test
  expect_error(check_num_between_0_and_1(1.00000001))
  expect_error(check_num_between_0_and_1(Inf))
  expect_error(check_num_between_0_and_1(100000))
})

test_that("check_num_between_0_and_1 gives errors for non-numeric inputs", {
  # Test
  expect_error(check_num_between_0_and_1("1"))
  expect_error(check_num_between_0_and_1(NA))
  expect_error(check_num_between_0_and_1(NULL))
  expect_error(check_num_between_0_and_1(matrix(1, ncol = 10, nrow = 10)))
  expect_error(check_num_between_0_and_1(rep(1, 10)))
})

# test check_if_dir_exists -----------------------------------------------------
test_that("check_if_dir_exists gives an error when dir doesn't exist", {
  # Set up
  temp_dir <- "/fake/directory/"

  # Test
  expect_error(check_if_dir_exists(temp_dir))
})

test_that("check_if_dir_exists doesn't give an error when dir does exist", {
  # Set up
  temp_dir <- "."

  # Test
  expect_error(check_if_dir_exists(temp_dir), NA)
})

# test check_if_permutation_num_valid ------------------------------------------
test_that("check_if_permutation_num_valid gives an error for a positive, non-integer perm = 1.5", {
  # Set up
  temp_perm <- 1.5

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error for a negative integer perm = -15", {
  # Set up
  temp_perm <- -15

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = 0", {
  # Set up
  temp_perm <- 0

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = '20' (a character)", {
  # Set up
  temp_perm <- '20'

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = NA", {
  # Set up
  temp_perm <- NA

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = NULL", {
  # Set up
  temp_perm <- NULL

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid doesn't give an error perm = 10000", {
  # Set up
  temp_perm <- 10000

  # Test
  expect_error(check_if_permutation_num_valid(temp_perm), NA)
})

# test check_is_string ---------------------------------------------------------
test_that("check_is_string gives an error when x = 100 (numeric)", {
  # Set up
  temp <- 100

  # Test
  expect_error(check_is_string(temp))
})

test_that("check_is_string gives an error when x = NA", {
  # Set up
  temp <- NA

  # Test
  expect_error(check_is_string(temp))
})

test_that("check_is_string gives an error when x = NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_is_string(temp))
})

test_that("check_is_string doesn't give an error when x = 'this is a string'", {
  # Set up
  temp <- 'this is a string'

  # Test
  expect_error(check_is_string(temp), NA)
})

# test check_if_vector ---------------------------------------------------------
test_that("check_if_vector doesn't give an error when x = c(1:10)", {
  # Set up
  temp <- c(1:10)

  # Test
  expect_error(check_if_vector(temp), NA)
})

test_that("check_if_vector doesn't give an error when x = letters[1:5]", {
  # Set up
  temp <- letters[1:5]

  # Test
  expect_error(check_if_vector(temp), NA)
})

test_that("check_if_vector gives an error when x = matrix(0, 10, 10)", {
  # Set up
  temp <- matrix(0, 10, 10)

  # Test
  expect_error(check_if_vector(temp))
})

test_that("check_if_vector gives an error when x = NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_if_vector(temp))
})

# test check_for_NA_and_inf ----------------------------------------------------
test_that("check_for_NA_and_inf gives an error when x is a dataframe", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing NA", {
  # Set up
  temp <- matrix(NA, 10, 10)

  # Test
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing -Inf", {
  # Set up
  temp <- matrix(-Inf, 10, 10)

  # Test
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing +Inf", {
  # Set up
  temp <- matrix(Inf, 10, 10)

  # Test
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf doesn't give an error when x is a matrix of zeroes", {
  # Set up
  temp <- matrix(0, 10, 10)

  # Test
  expect_error(check_for_NA_and_inf(temp), NA)
})

# test check_for_root_and_bootstrap
test_that("check_for_root_and_bootstrap doesn't give an error when x is rooted tree with node values of 100", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))

  # Tree
  expect_error(check_for_root_and_bootstrap(temp_tree), NA)
})

test_that("check_for_root_and_bootstrap gives an error when x is unrooted tree with node values of 100", {
  # Set up
  temp_tree <- rtree(20, rooted = FALSE)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node values of -100", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree))

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with too few node values", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree) - 1)

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with too many node values", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree) + 1)

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node values of NA", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(NA, Nnode(temp_tree))

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node.labels = NULL", {
  # Set up
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- NULL

  # Test
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when not given a tree object", {
  # Test
  expect_error(check_for_root_and_bootstrap("tree"))
  expect_error(check_for_root_and_bootstrap(10))
})

# test check_if_binary_vector --------------------------------------------------
test_that("check_if_binary_vector gives an error when x is rep(NA, 10)", {
  # Set up
  temp <- rep(NA, 10)

  # Test
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is letters[1:10]", {
  # Set up
  temp <- letters[1:10]

  # Test
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is dataframe(matrix(0, 10, 10))", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is matrix(0, 10, 10)", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector doesn't give an error when x is c(1, 0, 1, 0)", {
  # Set up
  temp <- c(1, 0, 1, 0)

  # Test
  expect_error(check_if_binary_vector(temp), NA)
})

test_that("check_if_binary_vector doesn't give an error when x is c(0, 0, 0, 0)", {
  # Set up
  temp <- c(0, 0, 0, 0)

  # Test
  expect_error(check_if_binary_vector(temp), NA)
})

test_that("check_if_binary_vector doesn't give an error when x is c(1)", {
  # Set up
  temp <- c(1)

  # Test
  expect_error(check_if_binary_vector(temp), NA)
})

# test check_if_binary_vector_numeric ------------------------------------------
test_that("check_if_binary_vector_numeric gives an error when x is rep(NA, 10)", {
  # Set up
  temp <- rep(NA, 10)

  # Test
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is letters[1:10]", {
  # Set up
  temp <- letters[1:10]

  # Test
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is dataframe(matrix(0, 10, 10))", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is matrix(0, 10, 10)", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(1, 0, 1, 0)", {
  # Set up
  temp <- c(1, 0, 1, 0)

  # Test
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(0, 0, 0, 0)", {
  # Set up
  temp <- c(0, 0, 0, 0)

  # Test
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(1)", {
  # Set up
  temp <- c(1)

  # Test
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

# test check_if_binary_matrix --------------------------------------------------
test_that("check_if_binary_matrix doesn't give an error when x is matrix(c(0, 1), 2, 1)", {
  # Set up
  temp <- matrix(c(0, 1), 2, 1)

  # Test
  expect_error(check_if_binary_matrix(temp), NA)
})

test_that("check_if_binary_matrix gives an error when x is matrix(NA, 2, 1)", {
  # Set up
  temp <- matrix(NA, 2, 1)

  # Test
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(Inf, 2, 1)", {
  # Set up
  temp <- matrix(Inf, 2, 1)

  # Test
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(c(1.5, 2.5), 2, 1)", {
  # Set up
  temp <- matrix(c(1.5, 2.5), 2, 1)

  # Test
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(3, 2, 1)", {
  # Set up
  temp <- matrix(3, 2, 1)

  # Test
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is as.data.frame(matrix(0, 2, 1))", {
  # Set up
  temp <- as.data.frame(matrix(0, 2, 1))

  # Test
  expect_error(check_if_binary_matrix(temp))
})

# test check_file_exists -------------------------------------------------------
test_that("check_file_exists gives an error when x is 'fake_file_name.txt'", {
  # Set up
  temp <- 'fake_file_name.txt'

  # Test
  expect_error(check_file_exists(temp))
})

test_that("check_file_exists doesn't give an error when x is 'test_validation.R'", {
  # Set up
  temp <- 'test_validation.R'

  # Test
  expect_error(check_file_exists(temp), NA)
})

# test check_rownames ----------------------------------------------------------
test_that("check_rownames doesn't give an error when tree$tip.label <- row.names(mat) <- letters[1:10]", {
  # Set up
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)
  temp_tree$tip.label <- row.names(temp_mat) <- letters[1:10]

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree), NA)
})

test_that("check_rownames gives an error when tree$tip.label <- letters[11:20],  row.names(mat) <- letters[1:10]", {
  # Set up
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)
  temp_tree$tip.label <- letters[11:20]
  row.names(temp_mat) <- letters[1:10]

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})

test_that("check_rownames gives an error when Ntip(tree) != nrow(mat)", {
  # Set up
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:8, nrow = 1, ncol = 8)

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))

  # Set up
  temp_tree <- rtree(2)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})


test_that("check_rownames gives an error when not given a tree object or not given matrix", {
  # Set up
  temp_tree <- rtree(10)
  fake_tree <- "tree"
  temp_mat <- matrix(1:10, nrow = 10, ncol = 1)
  row.names(temp_mat) <- temp_tree$tip.label
  fake_mat <- as.data.frame(temp_mat)

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree), NA)
  expect_error(check_rownames(mat = temp_mat, tr = fake_tree))
  expect_error(check_rownames(mat = fake_mat, tr = temp_tree))
})

test_that("check_rownames gives an error when given a matrix without rownames", {
  # Set up
  temp_tree <- rtree(10)
  temp_mat <- matrix(1:10, nrow = 10, ncol = 1)
  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))

  # Set up
  row.names(temp_mat) <- NULL

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})

test_that("check_rownames gives an error when given a tree without tip.labels", {
  # Set up
  temp_tree <- rtree(10)
  temp_mat <- matrix(1:10, nrow = 10, ncol = 1)
  temp_tree$tip.label <- NULL
  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})


test_that("check_rownames gives an error when tree$tip.label doesn't perfectly match matrix rownames", {
  # Set up
  set.seed(1)
  temp_tree <- rtree(10)
  temp_mat <- matrix(1:10, nrow = 10, ncol = 1)
  row.names(temp_mat) <- temp_tree$tip.label
  row.names(temp_mat)[1] <- "t6"
  row.names(temp_mat)[2] <- "t10"

  # Test
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})

# test check_is_number ---------------------------------------------------------
test_that("check_is_number doesn't give an error when x = 5", {
  # Set up
  temp <- 5

  # Test
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = 5.5", {
  # Set up
  temp <- 5.5

  # Test
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = pi", {
  # Set up
  temp <- pi

  # Test
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = -1/7", {
  # Set up
  temp <- -1/7

  # Test
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number gives an error when x = 'a'", {
  # Set up
  temp <- 'a'

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = matrix(0, 10, 10)", {
  # Set up
  temp <- matrix(0, 10, 10)

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = matrix(0, 10, 10)", {
  # Set up
  temp <- matrix(0, 10, 10)

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = data.frame(0, 10, 10)", {
  # Set up
  temp <- as.data.frame(matrix(0, 10, 10))

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = c(0, 1, 1)", {
  # Set up
  temp <- c(0, 1, 1)

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = NA", {
  # Set up
  temp <- NA

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = NULL", {
  # Set up
  temp <- NULL

  # Test
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = Inf", {
  # Set up
  temp <- Inf

  # Test
  expect_error(check_is_number(temp))
})

# test check_node_is_in_tree ---------------------------------------------------
test_that("check_node_is_in_tree doesn't give an error when node = 1", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = 1, tr = temp), NA)
})

test_that("check_node_is_in_tree doesn't give an error when node = Nnode(tree)", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = Nnode(temp), tr = temp), NA)
})

test_that("check_node_is_in_tree gives an error when node = 0", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = 0, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = -1", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = -1, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = Nnode + Ntip + 1", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = Nnode(temp) + Ntip(temp) + 1, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = 1.5", {
  # Set up
  temp <- rtree(10)
  temp$node.label <- c(1:Nnode(temp))

  # Test
  expect_error(check_node_is_in_tree(node_val = 1.5, tr = temp))
})

# test check_tree_is_valid -----------------------------------------------------
test_that("check_tree_is_valid returns true for randomly generated trees where Ntip is between 2 and 10", {
  # Test
  for (i in 2:10) {
    expect_error(check_tree_is_valid(rtree(i)), NA)
  }
})

test_that("check_tree_is_invalid throws an error when tree edge index is greater than Nedge(tree)", {
  # Set up
  invalid_tree <- rtree(10)
  for (j in 1:Nedge(invalid_tree)) {
    if (invalid_tree$edge[j, 2] == 1) {
      invalid_tree$edge[j, 2] <- Nedge(invalid_tree) + 1
      break
    }
  }

  # Test
  expect_error(check_tree_is_valid(invalid_tree))
})

# test check_convergence_possible ----------------------------------------------
test_that("check_convergence_possible gives an error 'discrete' and all values = 0", {
  # Set up
  disc_cont <- "discrete"
  temp_vec <- c(0,0,0,0,0,0,0,0)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible doesn't give an error 'discrete' and values= c(1,0,1,0,1,0,1,0)", {
  # Set up
  disc_cont <- "discrete"
  temp_vec <- c(1,0,1,0,1,0,1,0)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec), NA)
})


test_that("check_convergence_possible throws error given 'foobar' but valid values= c(1,0,1,0,1,0,1,0)", {
  # Set up
  disc_cont <- "foobar"
  temp_vec <- c(1,0,1,0,1,0,1,0)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible gives an error 'continuous' and all values = 0", {
  # Set up
  disc_cont <- "continuous"
  temp_vec <- c(0,0,0,0,0,0,0,0)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible gives an error 'discrete' and values= c(0, 0.1, 0.2)", {
  # Set up
  disc_cont <- "discrete"
  temp_vec <- c(0, 0.1, 0.2)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible gives an error 'discrete' and all values = 'a'", {
  # Set up
  disc_cont <- "discrete"
  temp_vec <- rep('a', 10)

  # Test
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

# test is_tip ------------------------------------------------------------------
test_that("is_tip returns TRUE when given a tree and the node = 1 (a tip)", {
  # Set up
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- 1

  # Test
  expect_true(is_tip(temp_node, temp_tree))
})

test_that("is_tip returns FALSE when given a tree and the node = Ntip(temp_tree) + 1 (not a tip)", {
  # Set up
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- Ntip(temp_tree) + 1

  # Test
  expect_false(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a tree and the node = 12.5", {
  # Set up
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- 12.5

  # Test
  expect_error(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a tree and the node = NA", {
  # Set up
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- NA

  # Test
  expect_error(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a matrix and a node", {
  # Set up
  temp_tree <- matrix(10, 10, 1)
  temp_node <- 5

  # Test
  expect_error(is_tip(temp_node, temp_tree))
})

# test check_if_g_mat_can_be_plotted -------------------------------------------
test_that("check_if_g_mat_can_be_plotted returns true for a binary matrix of 2x2", {
  # Set up
  temp_mat <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)

  # Test
  expect_error(check_if_g_mat_can_be_plotted(temp_mat), NA)
  expect_true(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns true for a binary matrix of 3x2 including NAs", {
  # Set up
  temp_mat <- matrix(c(0, 1, NA, 0, 1, NA), nrow = 3, ncol = 2)

  # Test
  expect_error(check_if_g_mat_can_be_plotted(temp_mat), NA)
  expect_true(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns an error for a non-binary matrix of 2x2", {
  # Set up
  temp_mat <- matrix(c(0.5, 1.5, 0, 1), nrow = 2, ncol = 2)

  # Test
  expect_error(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns FALSE for a binary matrix of all zeroes", {
  # Set up
  temp_mat <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)

  # Test
  expect_false(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns an FALSE for a binary matrix of 2x1", {
  # Set up
  temp_mat <- matrix(c(0, 1), nrow = 2, ncol = 1)

  # Test
  expect_false(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns FALSE for non-matrix, non-dataframe inputs", {
  # Test
  expect_false(check_if_g_mat_can_be_plotted(c(0, 1, 0, 1)))
  expect_false(check_if_g_mat_can_be_plotted(NA))
  expect_false(check_if_g_mat_can_be_plotted(NULL))
  expect_false(check_if_g_mat_can_be_plotted("temp"))
})

# test check_str_is_discrete_or_continuous -------------------------------------
test_that("check_str_is_discrete_or_continuous gives no error when given 'discrete'", {
  # Set up
  input <- 'discrete'

  # Test
  expect_error(check_str_is_discrete_or_continuous(input), NA)
})

test_that("check_str_is_discrete_or_continuous gives no error when given 'continuous'", {
  # Set up
  input <- 'continuous'

  # Test
  expect_error(check_str_is_discrete_or_continuous(input), NA)
})

test_that("check_str_is_discrete_or_continuous gives an error when given 'foobar'", {
  # Set up
  input <- 'foobar'

  # Test
  expect_error(check_str_is_discrete_or_continuous(input))
})

test_that("check_str_is_discrete_or_continuous gives an error when given NA", {
  # Set up
  input <- NA

  # Test
  expect_error(check_str_is_discrete_or_continuous(input))
})

test_that("check_str_is_discrete_or_continuous gives an error when given NULL", {
  # Set up
  input <- NULL

  # Test
  expect_error(check_str_is_discrete_or_continuous(input))
})

# test check_if_ancestral_reconstruction_method_compatible_with_ape ------------
test_that("check_if_ancestral_reconstruction_method_compatible_with_ape doesn't throw an error when given the acceptable inputs", {
  # Test
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape("ML"), NA)
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape("REML"), NA)
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape("pic"), NA)
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape("GLS"), NA)
})

test_that("check_if_ancestral_reconstruction_method_compatible_with_ape throws error when given the unacceptable inputs", {
  # Test
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape("other"))
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape(10))
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape(matrix(0, 10, 10)))
  expect_error(check_if_ancestral_reconstruction_method_compatible_with_ape(letters[1:10]))
})

# test check_equal -------------------------------------------------------------
test_that("check_equal doesn't throw errors when given two equal numbers", {
  # Test
  expect_error(check_equal(1, 1), NA)
  expect_error(check_equal(length(letters), 26), NA)
  expect_error(check_equal(20, 4 * 5), NA)
})

test_that("check_equal throws errors when given two not equal numbers", {
  # Test
  expect_error(check_equal(1, 0))
  expect_error(check_equal(length(letters), -26))
  expect_error(check_equal(20, 20.00000000001))
})

# test check_class -------------------------------------------------------------
test_that("check_class doesn't throw error when given an object and expected class", {
  # Test
  expect_error(check_class(1, "numeric"), NA)
  expect_error(check_class(c(1:10), "integer"), NA)
  expect_error(check_class("test", "character"), NA)
  expect_error(check_class(matrix(0), "matrix"), NA)
  expect_error(check_class(as.data.frame(matrix(0)), "data.frame"), NA)
  expect_error(check_class(rtree(2), "phylo"), NA)
})

test_that("check_class throws error when given an object and an incorrect class", {
  # Test
  expect_error(check_class(1, "integer"))
  expect_error(check_class(c(1:10), "character"))
  expect_error(check_class("test", "matrix"))
  expect_error(check_class(matrix(0), "data.frame"))
  expect_error(check_class(as.data.frame(matrix(0)), "phylo"))
  expect_error(check_class(rtree(2), "numeric"))
})
