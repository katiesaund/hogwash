library(microbeGWAS)

# TODO
# test check_input_format ------------------------------------------------------
#       This is a complex function with lots of inputs! This will take a bit or work to test!!!


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
  expect_error(check_if_vector(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing NA", {
  temp <- as.data.frame(matrix(NA, 10, 10))
  expect_error(check_if_vector(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing NULL", {
  temp <- as.data.frame(matrix(NULL, 10, 10))
  expect_error(check_if_vector(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing -Inf", {
  temp <- matrix(0, 10, 10)
  expect_error(check_if_vector(temp))
})

test_that("check_for_NA_and_inf gives an error when x is a matrix containing +Inf", {
  temp <- as.data.frame(matrix(Inf, 10, 10))
  expect_error(check_if_vector(temp))
})

test_that("check_for_NA_and_inf doesn't give an error when x is a matrix of zeroes", {
  temp <- as.data.frame(matrix(-Inf, 10, 10))
  expect_error(check_if_vector(temp))
})


#
# check_for_NA_and_inf <- function(mat){
#   # Function description -------------------------------------------------------
#   # Check that matrix contains no NAs and no +/- infinities.
#   #
#   # Input:
#   # mat. Matrix.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (class(mat) != "matrix"){
#     stop("Input should be a matrix.")
#   }
#   if (sum(is.na(mat)) > 0){
#     stop("Input matrices should not have any NA values.")
#   }
#   if (sum(mat == -Inf) > 0){
#     stop("Inpute matrices should not have any -Inf values.")
#   }
#   if (sum(mat == Inf) > 0){
#     stop("Inpute matrices should not have any -Inf values.")
#   }
# } # end check_for_NA_and_inf()
#
# check_for_root_and_boostrap <- function(tr){
#   # Function description -------------------------------------------------------
#   # Check that phylogenetic tree is rooted and contains bootstrap values in the node labels.
#   #
#   # Input:
#   # tr. Phylo.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (class(tr) != "phylo"){
#     stop("Tree must be phylo object")
#   }
#   if (!is.rooted(tr)){
#     stop("Tree must be rooted")
#   }
#   if (is.null(tr$node.label)){
#     stop("Tree must have bootstrap values in the nodes")
#   }
# } # end check_for_root_and_boostrap()
#
# check_if_binary_vector <- function(vec){
#   # Function description -------------------------------------------------------
#   # Check that the matrix only contains values 1 or 0.
#   #
#   # Input:
#   # vec. Vector.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (sum(!(vec %in% c(0, 1))) > 0 | class(vec) != "integer"){
#     stop("Vector should be only 1s and 0s")
#   }
# } # end check_if_binary_vector()
#
# check_if_binary_vector_numeric <- function(vec){
#   # Function description -------------------------------------------------------
#   # Check that the matrix only contains values 1 or 0.
#   #
#   # Input:
#   # vec. Vector.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (sum(!(vec %in% c(0, 1))) > 0 | class(vec) != "numeric"){
#     stop("Vector should be only 1s and 0s")
#   }
# } # end check_if_binary_vector_numeric()
#
#
# check_if_binary_matrix <- function(mat){
#   # Function description -------------------------------------------------------
#   # Check that the matrix only contains values 1 or 0.
#   #
#   # Input:
#   # mat. Matrix.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (sum(!(mat %in% c(0, 1))) > 0 | class(mat) != "matrix"){
#     stop("Genotype matrix should be only 1s and 0s")
#   }
# } # end check_if_binary_matrix()
#
# check_file_exists <- function(file_name){
#   # Function description -------------------------------------------------------
#   # Check that the file exists.
#   #
#   # Input:
#   # file_name. Character.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (!file.exists(file_name)){
#     stop("File does not exist")
#   }
# } # end check_file_exists()
#
# check_rownames <- function(mat, tr){
#   # Function description -------------------------------------------------------
#   # Check that phylogenetic tree tip labels are identical to the matrix row.names.
#   #
#   # Input:
#   # mat. Matrix.
#   # tr. Phylo.
#   #
#   # Output:
#   # None.
#   #
#   # Check input ----------------------------------------------------------------
#   if (class(mat) != "matrix" | class(tr) != "phylo"){
#     stop("Inputs are incorrectly formatted.")
#   }
#
#   # Function -------------------------------------------------------------------
#   if (sum(row.names(mat) != tr$tip.label) != 0){
#     stop("Matrix must be formatted with samples in matrix in the same order as tree$tip.label.")
#   }
# } # end check_rownames()
#
# check_is_number <- function(num){
#   # Function description -------------------------------------------------------
#   # Check that input is some type of number.
#   #
#   # Input:
#   # num. Number. Could be numeric, double, or integer.
#   #
#   # Output:
#   # None.
#   #
#   # Check input & function -----------------------------------------------------
#   if (!is.numeric(num)){
#     if (!is.integer(num)){
#       if (!is.double(num)){
#         stop("Must be a number")
#       }
#     }
#   }
# } # end check_is_number()
#
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
