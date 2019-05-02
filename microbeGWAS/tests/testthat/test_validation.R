library(microbeGWAS)

# TODO
# test check_input_format ------------------------------------------------------
#       This is a complex function with lots of inputs! This will take a bit or work to test!!!

# TODO check_is_number -- Should I also add checks to remove Inf? Inf is a numeric...but I can't imagine that would be good for GWAS...
# TODO check_node_is_in_tree -- what is the number of nodes given number of tips in a bifurcating tree? use this value to update the failing test

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

test_that("check_node_is_in_tree gives an error when node = Nnode + 1", {
  temp = rtree(10)
  temp$node.label <- c(1:Nnode(temp))
  expect_error(check_node_is_in_tree(node_val = Nnode(temp) + 1, tr = temp))
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

# check_node_is_in_tree <- function(node_val, tr){
#   # Function description ------------------------------------------------------#
#   # Test if a node value is plausibly contained within the tree.
#   #
#   # Inputs:
#   # node_val: Integer. Index of node.
#   # tr: phylogenetic tree.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function ----------------------------------------------------#
#   check_for_root_and_boostrap(tr)
#   check_is_number(node_val)
#
#   if (node_val > Nnode(tr) + Ntip(tr)){
#     stop("Node number is too high; not found in tree.")
#   }
#   if (node_val < 1 | !is.integer(node_val)){
#     stop("Node number must be positive integer")
#   }
# } # end check_node_is_in_tree()
#
# test_that("check_tree_is_valid returns true for randomly generated trees where Ntip is between 2 and 10", {
#   for (i in 2:10){
#     expect_true(check_tree_is_valid(rtree(i)))
#   }
#
#
# })
#
# test_that("check_tree_is_invalid throws an error when tree edge index is greater than Nedge(tree)", {
#   invalid_tree <- rtree(10)
#   for (j in 1:Nedge(invalid_tree)){
#     if (invalid_tree$edge[j, 2] == 1){
#       invalid_tree$edge[j, 2] <- Nedge(invalid_tree) + 1
#       break
#     }
#   }
#   expect_that(check_tree_is_valid(invalid_tree), throws_error())
# })
#
# check_tree_is_valid <- function(tr){
#   #` A valid tree has N nodes, with n_tips nodes being tip nodes, numbered 1 through n_tips`
#
#   num_edges_for_node <- table(tr$edge)
#
#   for (i in 1:Ntip(tr)){
#     if (num_edges_for_node[i] != 1){
#       stop(paste("Tip node", i, "has", num_edges_for_node[i], "edges. Should have 1 edge"))
#     }
#   }
#   for (i in (Ntip(tr) + 1):(Nnode(tr) + Ntip(tr))){
#     if (num_edges_for_node[i] != 2 && num_edges_for_node[i] != 3){
#       stop(paste("Internal node", i, "has", num_edges_for_node[i], "edges. Should have 2(root) or 3 edge"))
#     }
#   }
#   return(TRUE)
# }
#
# # check_transitions <- function(transition_vector, tr){
# #   for (i in 1:length(transition_vector)){
# #     parent <- tr$edge[i, 1]
# #     child <- tr$edge[i, 2]
# #
# #   }
# #
# # } # end check_transitions()
#
# check_convergence_possible <- function(vec, discrete_or_continuous){
#   convergence_not_possible <- FALSE
#   if (discrete_or_continuous == "discrete"){
#     check_if_binary_vector(vec)
#     if (sum(vec) >= (length(vec)-1) | sum(vec) <= 1){
#       convergence_not_possible <- TRUE
#     }
#
#   } else if (discrete_or_continuous == "continuous"){
#     if (length(unique(vec)) == 1){
#       convergence_not_possible <- TRUE
#     }
#   }
#   if (convergence_not_possible){
#     stop("Convergence is not possible for this phenotype")
#   }
# } # end check_convergence_possible()
#
#
#
#
# context("Complex functions") #-------------------------------------------------#
# test_that("ancestral_reconstruction_by_ML with discrete input produce ancestral reconstruction with a value for each tip and node.", {
#   tree <- rtree(9, rooted = TRUE)
#   tree$node.label <- rep(100, Nnode(tree))
#   num_col <- 9
#   num_cells <- num_col * Ntip(tree)
#   test_mat <- matrix(rep(c(1, 0), num_cells), nrow = Ntip(tree), ncol = num_col)
#   check_if_binary_matrix(test_mat)
#   dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "discrete")
#   dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "discrete")
#   expected_length <- Ntip(tree) + Nnode(tree)
#   expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
#   expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
# })
#
# test_that("ancestral_reconstruction_by_ML with continuous input produce ancestral reconstruction with a value for each tip and node.", {
#   tree <- rtree(9, rooted = TRUE)
#   tree$node.label <- rep(100, Nnode(tree))
#   num_col <- 9
#   num_cells <- num_col * Ntip(tree)
#   test_mat <- matrix(rnorm(num_cells, mean = 0, sd = 10), nrow = Ntip(tree), ncol = num_col)
#   dummy_pheno <- ancestral_reconstruction_by_ML(tree, test_mat, 1, "continuous")
#   dummy_geno <-  ancestral_reconstruction_by_ML(tree, test_mat, 2, "continuous")
#   expected_length <- Ntip(tree) + Nnode(tree)
#   expect_identical(length(dummy_pheno$tip_and_node_recon), expected_length)
#   expect_identical(length(dummy_geno$tip_and_node_recon), expected_length)
# })
