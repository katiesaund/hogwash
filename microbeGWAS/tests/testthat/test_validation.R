library(microbeGWAS)
context("Validation functions") #----------------------------------------------#

# test check_dimensions --------------------------------------------------------
test_that("check_dimensions gives an error for when given a dataframe instead of a matrix", {
  temp_mat <- as.data.frame(matrix(1:100, nrow = 10, ncol = 10))
  expect_error(check_dimensions(temp_mat, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given a vector instead of a matrix", {
  temp_vector <- c(1:100)
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given min_rows is too large", {
  num_row <- 10
  temp_mat <- matrix(rep("A", 100), nrow = num_row, ncol = 10)
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = 10,
                                min_rows = num_row + 99,
                                exact_cols = NULL,
                                min_cols = 1))
})

test_that("check_dimensions gives an error for when given a min_cols it too large", {
  num_col <- 10
  temp_mat <- matrix(rep("A", 100), nrow = 10, ncol = num_col)
  expect_error(check_dimensions(mat = temp_mat,
                                exact_rows = NULL,
                                min_rows = 1,
                                exact_cols = NULL,
                                min_cols = num_col + 2))
})

test_that("check_dimensions gives an error for when given exact_cols that is smaller or larger than ncol()", {
  num_col <- 10
  temp_mat <- matrix(rep("A", 100), nrow = 2, ncol = num_col)
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
  num_row <- 10
  temp_mat <- matrix(rep("A", 100), nrow = num_row, ncol = 1)
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
  temp_vector <- NULL
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions gives an error for when given NA instead of a matrix", {
  temp_vector <- NULL
  expect_error(check_dimensions(temp_vector, 10, 1, 10, 1))
})

test_that("check_dimensions doesn't give an error for when given a valid matrix", {
  temp_matrix <- matrix(0, nrow = 10, ncol = 5)
  expect_error(check_dimensions(temp_vector, 10, 1, 5, 1))
})

# test check_if_alpha_valid ----------------------------------------------------
test_that("check_if_alpha_valid gives an error alpha = 5", {
  temp_alpha <- 5
  expect_error(check_if_alpha_valid(temp_alpha))
})

test_that("check_if_alpha_valid gives an error alpha = -0.1", {
  temp_alpha <- -0.1
  expect_error(check_if_alpha_valid(temp_alpha))
})

test_that("check_if_alpha_valid gives an error alpha = 'A'", {
  temp_alpha <- 'A'
  expect_error(check_if_alpha_valid(temp_alpha))
})

test_that("check_if_alpha_valid gives an error when alpha is a matrix", {
  temp_alpha <- matrix(0.5)
  expect_error(check_if_alpha_valid(temp_alpha))
})


# test check_if_dir_exists -----------------------------------------------------
test_that("check_if_dir_exists gives an error when dir doesn't exist", {
  temp_dir <- "/fake/directory/"
  expect_error(check_if_dir_exists(temp_dir))
})

test_that("check_if_dir_exists doesn't give an error when dir does exist", {
  temp_dir <- "."
  expect_error(check_if_dir_exists(temp_dir), NA)
})

# test check_if_permutation_num_valid ------------------------------------------
test_that("check_if_permutation_num_valid gives an error whem perm = 1.5", {
  temp_perm <- 1.5
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = -1.5", {
  temp_perm <- -1.5
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = 0", {
  temp_perm <- 0
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = '20' (a char)", {
  temp_perm <- '20'
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = NA", {
  temp_perm <- NA
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid gives an error whem perm = NULL", {
  temp_perm <- NULL
  expect_error(check_if_permutation_num_valid(temp_perm))
})

test_that("check_if_permutation_num_valid doesn't give an error perm = 10000", {
  temp_perm <- 10000
  expect_error(check_if_permutation_num_valid(temp_perm), NA)
})

# test check_is_string ---------------------------------------------------------
test_that("check_is_string gives an error when x = 100 (numeric)", {
  temp <- 100
  expect_error(check_is_string(temp))
})

test_that("check_is_string gives an error when x = NA", {
  temp <- NA
  expect_error(check_is_string(temp))
})

test_that("check_is_string gives an error when x = NULL", {
  temp <- NULL
  expect_error(check_is_string(temp))
})

test_that("check_is_string doesn't give an error when x = 'this is a string'", {
  temp <- 'this is a string'
  expect_error(check_is_string(temp), NA)
})


test_that("check_is_string doesn't give an error when x = 'this is a string'", {
  temp <- 'this is a string'
  expect_error(check_is_string(temp), NA)
})

# test check_if_vector ---------------------------------------------------------
test_that("check_if_vector doesn't give an error when x = c(1:10)", {
  temp <- c(1:10)
  expect_error(check_if_vector(temp), NA)
})

test_that("check_if_vector doesn't give an error when x = letters[1:5]", {
  temp <- letters[1:5]
  expect_error(check_if_vector(temp), NA)
})

test_that("check_if_vector gives an error when x = matrix(0, 10, 10)", {
  temp <- matrix(0, 10, 10)
  expect_error(check_if_vector(temp))
})

test_that("check_if_vector gives an error when x = NULL", {
  temp <- NULL
  expect_error(check_if_vector(temp))
})

# test check_for_NA_and_inf ----------------------------------------------------
test_that("check_for_NA_and_inf gives an error when x is a dataframe", {
  temp <- as.data.frame(matrix(0, 10, 10))
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing NA", {
  temp <- matrix(NA, 10, 10)
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is NULL", {
  temp <- NULL
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing -Inf", {
  temp <- matrix(-Inf, 10, 10)
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing +Inf", {
  temp <- matrix(Inf, 10, 10)
  expect_error(check_for_NA_and_inf(temp))
})

test_that("check_for_NA_and_inf doesn't give an error when x is a matrix of zeroes", {
  temp <- matrix(0, 10, 10)
  expect_error(check_for_NA_and_inf(temp), NA)
})

# test check_for_root_and_bootstrap
test_that("check_for_root_and_bootstrap doesn't give an error when x is rooted tree with node values of 100", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  expect_error(check_for_root_and_bootstrap(temp_tree), NA)
})

test_that("check_for_root_and_bootstrap gives an error when x is unrooted tree with node values of 100", {
  temp_tree <- rtree(20, rooted = FALSE)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node values of -100", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree))
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with too few node values", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree) - 1)
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with too many node values", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(-100, Nnode(temp_tree) + 1)
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node values of NA", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- rep(NA, Nnode(temp_tree))
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

test_that("check_for_root_and_bootstrap gives an error when x is rooted tree with node.labels = NULL", {
  temp_tree <- rtree(20, rooted = TRUE)
  temp_tree$node.labels <- NULL
  expect_error(check_for_root_and_bootstrap(temp_tree))
})

# test check_if_binary_vector --------------------------------------------------
test_that("check_if_binary_vector gives an error when x is rep(NA, 10)", {
  temp <- rep(NA, 10)
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is NULL", {
  temp <- NULL
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is letters[1:10]", {
  temp <- letters[1:10]
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is dataframe(matrix(0, 10, 10))", {
  temp <- as.data.frame(matrix(0, 10, 10))
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector gives an error when x is matrix(0, 10, 10)", {
  temp <- as.data.frame(matrix(0, 10, 10))
  expect_error(check_if_binary_vector(temp))
})

test_that("check_if_binary_vector doesn't give an error when x is c(1, 0, 1, 0)", {
  temp <- c(1, 0, 1, 0)
  expect_error(check_if_binary_vector(temp), NA)
})

test_that("check_if_binary_vector doesn't give an error when x is c(0, 0, 0, 0)", {
  temp <- c(0, 0, 0, 0)
  expect_error(check_if_binary_vector(temp), NA)
})

test_that("check_if_binary_vector doesn't give an error when x is c(1)", {
  temp <- c(1)
  expect_error(check_if_binary_vector(temp), NA)
})

# test check_if_binary_vector_numeric ------------------------------------------
test_that("check_if_binary_vector_numeric gives an error when x is rep(NA, 10)", {
  temp <- rep(NA, 10)
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is NULL", {
  temp <- NULL
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is letters[1:10]", {
  temp <- letters[1:10]
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is dataframe(matrix(0, 10, 10))", {
  temp <- as.data.frame(matrix(0, 10, 10))
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric gives an error when x is matrix(0, 10, 10)", {
  temp <- as.data.frame(matrix(0, 10, 10))
  expect_error(check_if_binary_vector_numeric(temp))
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(1, 0, 1, 0)", {
  temp <- c(1, 0, 1, 0)
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(0, 0, 0, 0)", {
  temp <- c(0, 0, 0, 0)
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

test_that("check_if_binary_vector_numeric doesn't give an error when x is c(1)", {
  temp <- c(1)
  expect_error(check_if_binary_vector_numeric(temp), NA)
})

# test check_if_binary_matrix --------------------------------------------------
test_that("check_if_binary_matrix doesn't give an error when x is matrix(c(0, 1), 2, 1)", {
  temp <- matrix(c(0, 1), 2, 1)
  expect_error(check_if_binary_matrix(temp), NA)
})

test_that("check_if_binary_matrix gives an error when x is matrix(NA, 2, 1)", {
  temp <- matrix(NA, 2, 1)
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(Inf, 2, 1)", {
  temp <- matrix(Inf, 2, 1)
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is NULL", {
  temp <- NULL
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(c(1.5, 2.5), 2, 1)", {
  temp <- matrix(c(1.5, 2.5), 2, 1)
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is matrix(3, 2, 1)", {
  temp <- matrix(3, 2, 1)
  expect_error(check_if_binary_matrix(temp))
})

test_that("check_if_binary_matrix gives an error when x is as.data.frame(matrix(0, 2, 1))", {
  temp <- as.data.frame(matrix(0, 2, 1))
  expect_error(check_if_binary_matrix(temp))
})

# test check_file_exists -------------------------------------------------------
test_that("check_file_exists gives an error when x is 'fake_file_name.txt'", {
  temp <- 'fake_file_name.txt'
  expect_error(check_file_exists(temp))
})

test_that("check_file_exists doesn't give an error when x is 'test_validation.R'", {
  temp <- 'test_validation.R'
  expect_error(check_file_exists(temp), NA)
})

# test check_rownames ----------------------------------------------------------
test_that("check_rownames doesn't give an error when tree$tip.label <- row.names(mat) <- letters[1:10]", {
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)
  temp_tree$tip.label <- row.names(temp_mat) <- letters[1:10]
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree), NA)
})

test_that("check_rownames gives an error when tree$tip.label <- letters[11:20],  row.names(mat) <- letters[1:10]", {
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)
  temp_tree$tip.label <- letters[11:20]
  row.names(temp_mat) <- letters[1:10]
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))
})

test_that("check_rownames gives an error when Ntip(tree) != nrow(mat)", {
  temp_tree <- rtree(10)
  temp_mat  <- matrix(1:8, nrow = 1, ncol = 8)
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))

  temp_tree <- rtree(2)
  temp_mat  <- matrix(1:80, nrow = 10, ncol = 8)
  expect_error(check_rownames(mat = temp_mat, tr = temp_tree))

})

# test check_is_number ---------------------------------------------------------
test_that("check_is_number doesn't give an error when x = 5", {
  temp = 5
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = 5.5", {
  temp = 5.5
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = pi", {
  temp = pi
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number doesn't give an error when x = -1/7", {
  temp = -1/7
  expect_error(check_is_number(temp), NA)
})

test_that("check_is_number gives an error when x = 'a'", {
  temp = 'a'
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = matrix(0, 10, 10)", {
  temp = matrix(0, 10, 10)
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = matrix(0, 1, 1)", {
  temp = matrix(0, 10, 10)
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = c(0, 1, 1)", {
  temp = c(0, 1, 1)
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = NA", {
  temp = NA
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = NULL", {
  temp = NULL
  expect_error(check_is_number(temp))
})

test_that("check_is_number gives an error when x = Inf", {
  temp = Inf
  expect_error(check_is_number(temp))
})

# test check_node_is_in_tree ---------------------------------------------------
test_that("check_node_is_in_tree doesn't give an error when node = 1", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = 1, tr = temp), NA)
})

test_that("check_node_is_in_tree doesn't give an error when node = Nnode(tree)", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = Nnode(temp), tr = temp), NA)
})

test_that("check_node_is_in_tree gives an error when node = 0", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = 0, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = -1", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = -1, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = Nnode + Ntip + 1", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = Nnode(temp) + Ntip(temp) + 1, tr = temp))
})

test_that("check_node_is_in_tree gives an error when node = 1.5", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = 1.5, tr = temp))
})


# test check_tree_is_valid -----------------------------------------------------
test_that("check_tree_is_valid returns true for randomly generated trees where Ntip is between 2 and 10", {
  for (i in 2:10){
    expect_error(check_tree_is_valid(rtree(i)), NA)
  }
})

test_that("check_tree_is_invalid throws an error when tree edge index is greater than Nedge(tree)", {
  invalid_tree <- rtree(10)
  for (j in 1:Nedge(invalid_tree)){
    if (invalid_tree$edge[j, 2] == 1){
      invalid_tree$edge[j, 2] <- Nedge(invalid_tree) + 1
      break
    }
  }
  expect_error(check_tree_is_valid(invalid_tree))
})

# test check_convergence_possible ----------------------------------------------
test_that("check_convergence_possible gives an error 'discrete' and all values = 0", {
  disc_cont <- "discrete"
  temp_vec <- c(0,0,0,0,0,0,0,0)
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible doesn't give an error 'discrete' and values= c(1,0,1,0,1,0,1,0)", {
  disc_cont <- "discrete"
  temp_vec <- c(1,0,1,0,1,0,1,0)
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec), NA)
})

test_that("check_convergence_possible gives an error 'continuous' and all values = 0", {
  disc_cont <- "continuous"
  temp_vec <- c(0,0,0,0,0,0,0,0)
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

test_that("check_convergence_possible gives an error 'discrete' and values= c(0, 0.1, 0.2)", {
  disc_cont <- "discrete"
  temp_vec <- c(0, 0.1, 0.2)
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})


test_that("check_convergence_possible gives an error 'discrete' and all values = 'a'", {
  disc_cont <- "discrete"
  temp_vec <- rep('a', 10)
  expect_error(check_convergence_possible(discrete_or_continuous = disc_cont, vec = temp_vec))
})

# test is_tip ------------------------------------------------------------------
test_that("is_tip returns TRUE when given a tree and the node = 1 (a tip)", {
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- 1
  expect_true(is_tip(temp_node, temp_tree))
})

test_that("is_tip returns FALSE when given a tree and the node = Ntip(temp_tree) + 1 (not a tip)", {
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- Ntip(temp_tree) + 1
  expect_false(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a tree and the node = 12.5", {
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- 12.5
  expect_error(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a tree and the node = NA", {
  temp_tree <- rtree(10)
  temp_tree$node.labels <- rep(100, Nnode(temp_tree))
  temp_node <- NA
  expect_error(is_tip(temp_node, temp_tree))
})

test_that("is_tip gives an error when given a matrix and a node", {
  temp_tree <- matrix(10, 10, 1)
  temp_node <- 5
  expect_error(is_tip(temp_node, temp_tree))
})

# test check_if_g_mat_can_be_plotted -------------------------------------------
test_that("check_if_g_mat_can_be_plotted returns true for a binary matrix of 2x2", {
  temp_mat <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  expect_true(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns an error for a non-binary matrix of 2x2", {
  temp_mat <- matrix(c(0.5, 1.5, 0, 1), nrow = 2, ncol = 2)
  expect_error(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns an error for a binary dataframe of 2x2", {
  temp_mat <- as.data.frame(matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2))
  expect_error(check_if_g_mat_can_be_plotted(temp_mat))
})

test_that("check_if_g_mat_can_be_plotted returns FALSE for a binary matrix of all zeroes", {
  temp_mat <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  expect_false(check_if_g_mat_can_be_plotted(temp_mat))
})

# check_if_g_mat_can_be_plotted <- function(geno_matrix){
#   # Function description -------------------------------------------------------
#   # TODO
#   # In order to make a heatmap there need to be 1) at least two columns, 2)
#   # two different values within the matrix (0 and 1).
#   # There can be NAs in the matrix.
#   #
#   # Inputs:
#   # geno_matrix. Matrix. 1, 0, and/or NA.
#   #
#   # Outputs:
#   # plot_logical. Logical. TRUE or FALSE.
#   #
#   # Check input ----------------------------------------------------------------
#   check_dimensions(geno_matrix, min_rows = 1, min_cols = 2)
#
#   if (sum(as.vector(geno_matrix)[!is.na(as.vector(geno_matrix))] %% 1 == 0) != 0){
#     stop("Joint genotype matrix + phenotype must contain only 1, 0, or NA. (For discrete heatmap plot).")
#   }
#
#   # Function -------------------------------------------------------------------
#   ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
#   zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
#   nas <- sum(is.na(geno_matrix)) > 1
#
#   plot_logical <- FALSE #
#   if (ones == 1 && zeros == 1 && nas == 0) {
#     plot_logical <- TRUE
#   }
#   if (ones + zeros + nas == 3) {
#     plot_logical <- TRUE
#   }
#
#   # Return output --------------------------------------------------------------
#   if(!is.logical(plot_logical)){stop("Output must be a logical")}
#   return(plot_logical)
# }
