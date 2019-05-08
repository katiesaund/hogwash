# Library of all functions that "check" or "assert" that something is true.
# If any functions get added to this file, be sure to add test() to test_validation.R

# TODO Add function description for check_if_g_mat_can_be_plotted

check_dimensions <- function(mat, exact_rows = NULL, min_rows, exact_cols = NULL, min_cols){
  # Function description -------------------------------------------------------
  # Check that the matrix is of the specified dimensions.
  #
  # Input:
  # mat:        Matrix.
  # exact_rows. Numeric. Describes expected number of rows in matrix. Default is NULL.
  # min_rows.   Numeric. Describes minimum number of rows in matrix. Must be specified.
  # exact_cols. Numeric. Describes expected number of columns in matrix. Default NULL.
  # min_cols.   Numeric. Describes minimum number of columns in matrix. Must be specified.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  if (!is.null(exact_rows)){check_is_number(exact_rows)}
  check_is_number(min_rows)
  if (!is.null(exact_cols)){check_is_number(exact_cols)}
  check_is_number(min_cols)
  if (class(mat) != "matrix"){stop("input must be a matrix")}

  # Function -------------------------------------------------------------------
  if (nrow(mat) < min_rows){
    stop("matrix has too few rows")
  }
  if (ncol(mat) < min_cols){
    stop("matrix has too few columns")
  }
  if (!(is.null(exact_rows))){
    if (nrow(mat) != exact_rows){
      stop("matrix has wrong number of rows")
    }
  }
  if (!(is.null(exact_cols))){
    if (ncol(mat) != exact_cols){
      stop("matrix has wrong number of columns")
    }
  }
} # end check_dimensions()

check_if_alpha_valid <- function(a){
  # Function description -------------------------------------------------------
  # Check that the alpha (threshold for significance) is within a valid range (0 < alpha < 1).
  #
  # Input:
  # a. Numeric. Value of alpha.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_is_number(a)

  # Function -------------------------------------------------------------------
  if (a > 1 | a < 0){
    stop("Provide a valid alpha.")
  }
} # end check_if_alpha_valid()

check_if_dir_exists <- function(dir){
  # Function description -------------------------------------------------------
  # Check that output directory exists so data can be saved in it.
  #
  # Input:
  # dir. Character. Path to output directory.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_is_string(dir)

  # Function -------------------------------------------------------------------
  if (!dir.exists(dir)){
    stop("The output directory indicated does not exist.")
  }
} # end check_if_dir_exists()

check_if_permutation_num_valid <- function(perm){
  # Function description -------------------------------------------------------
  # Check that the permutation number inidicated is valid (1 <= perm).
  #
  # Input:
  # perm. Number. Times to shuffle the data on the tree to create a null distribution for the permutation test.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_is_number(perm)

  # Function -------------------------------------------------------------------
  if (perm < 1 ||perm %% 1 != 0){
    stop("The permutation number should be a positive integer indicating the number of null distributions to create.")
  }
} # end check_if_permutation_num_valid()

check_is_string <- function(char){
  # Function description -------------------------------------------------------
  # Check that the input is a character string.
  #
  # Input:
  # char. Character.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (!is.character(char)){
    stop("Object must be a character string.")
  }
} # end check_is_string()

check_if_vector <- function(vector){
  # Function description -------------------------------------------------------
  # Check that input is a vector.
  #
  # Input:
  # vector. Vector.

  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (!is.vector(vector)){
    stop("Input must be a vector")
  }
} # end check_if_vector()

check_for_NA_and_inf <- function(mat){
  # Function description -------------------------------------------------------
  # Check that matrix contains no NAs and no +/- infinities.
  #
  # Input:
  # mat. Matrix.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (class(mat) != "matrix"){
    stop("Input should be a matrix.")
  }
  if (sum(is.na(mat)) > 0){
    stop("Input matrices should not have any NA values.")
  }
  if (sum(mat == -Inf) > 0){
    stop("Inpute matrices should not have any -Inf values.")
  }
  if (sum(mat == Inf) > 0){
    stop("Inpute matrices should not have any -Inf values.")
  }
} # end check_for_NA_and_inf()

check_for_root_and_bootstrap <- function(tr){
  # Function description -------------------------------------------------------
  # Check that phylogenetic tree is rooted and contains bootstrap values in the node labels.
  #
  # Input:
  # tr. Phylo.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (class(tr) != "phylo"){
    stop("Tree must be phylo object")
  }
  if (!is.rooted(tr)){
    stop("Tree must be rooted")
  }
  if (is.null(tr$node.label)){
    stop("Tree must have bootstrap values in the nodes")
  }
  if (length(tr$node.label) != Nnode(tr)){
    stop("Tree must have bootstrap values for each node")

  }
  for (i in 1:length(tr$node.label)){
    check_is_number(tr$node.label[i]) # Node lables must be bootstrap values
    if (tr$node.label[i] < 0){
      stop("Tree bootstrap values must be >= 0")
    }
  }
} # end check_for_root_and_bootstrap()

check_if_binary_vector <- function(vec){
  # Function description -------------------------------------------------------
  # Check that the matrix only contains values 1 or 0.
  #
  # Input:
  # vec. Vector.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (sum(!(vec %in% c(0, 1))) > 0){
    stop("Vector should be only 1s and 0s")
  }
  if (class(vec) != "integer"){
    if (class(vec) != "numeric"){
      stop("Vector should be only 1s and 0s")
    }
  }
} # end check_if_binary_vector()

check_if_binary_vector_numeric <- function(vec){
  # Function description -------------------------------------------------------
  # Check that the matrix only contains values 1 or 0.
  #
  # Input:
  # vec. Vector.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (sum(!(vec %in% c(0, 1))) > 0 | class(vec) != "numeric"){
    stop("Vector should be only 1s and 0s")
  }
} # end check_if_binary_vector_numeric()


check_if_binary_matrix <- function(mat){
  # Function description -------------------------------------------------------
  # Check that the matrix only contains values 1 or 0.
  #
  # Input:
  # mat. Matrix.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (sum(!(mat %in% c(0, 1))) > 0 | class(mat) != "matrix"){
    stop("Genotype matrix should be only 1s and 0s")
  }
  if (is.null(dim(mat))){
    stop("Must be a matrix")
  }
  if (nrow(mat) == 0 | ncol(mat) == 0){
    stop("matrix must have columns and rows.")
  }
} # end check_if_binary_matrix()

check_file_exists <- function(file_name){
  # Function description -------------------------------------------------------
  # Check that the file exists.
  #
  # Input:
  # file_name. Character.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (!file.exists(file_name)){
    stop("File does not exist")
  }
} # end check_file_exists()

check_rownames <- function(mat, tr){
  # Function description -------------------------------------------------------
  # Check that phylogenetic tree tip labels are identical to the matrix row.names.
  #
  # Input:
  # mat. Matrix.
  # tr. Phylo.
  #
  # Output:
  # None.
  #
  # Check input ----------------------------------------------------------------
  if (class(mat) != "matrix" | class(tr) != "phylo"){
    stop("Inputs are incorrectly formatted.")
  }

  # Function -------------------------------------------------------------------
  if(is.null(row.names(mat))){
    stop("Matrix must have row.names()")
  }
  if(is.null(tr$tip.label)){
    stop("Tree must have tip labels")
  }
  if (sum(row.names(mat) != tr$tip.label) != 0){
    stop("Matrix must be formatted with samples in matrix in the same order as tree$tip.label.")
  }
  if (sum(row.names(mat) == tr$tip.label) != Ntip(tr)){
    stop("Matrix must be formatted with samples in matrix in the same order as tree$tip.label.")
  }
} # end check_rownames()

check_is_number <- function(num){
  # Function description -------------------------------------------------------
  # Check that input is some type of number.
  # Check that the input is a single number.
  #
  # Input:
  # num. Number. Could be numeric, double, or integer.
  #
  # Output:
  # None.
  #
  # Check input & function -----------------------------------------------------
  if (!is.numeric(num)){
    if (!is.integer(num)){
      if (!is.double(num)){
        stop("Must be a number")
      }
    }
  }
  if (is.data.frame(num)){
    stop("Number can't be a dataframe")
  }
  if (is.matrix(num)){
    stop("Number can't be a matrix")
  }
  if (length(num) != 1){
    stop("Must be a single number")
  }
  if (num == Inf | num == -Inf){
    stop("This can't handle infinity.")
  }
} # end check_is_number()

check_node_is_in_tree <- function(node_val, tr){
  # Function description ------------------------------------------------------#
  # Test if a node value is plausibly contained within the tree.
  #
  # Inputs:
  # node_val: Integer. Index of node.
  # tr: phylogenetic tree.
  #
  # Output:
  # None.
  #
  # Check input & function ----------------------------------------------------#
  check_for_root_and_bootstrap(tr)
  check_is_number(node_val)

  if (node_val > Nnode(tr) + Ntip(tr)){
    stop("Node number is too high; not found in tree.")
  }
  if (node_val < 1 | node_val %% 1 != 0){
    stop("Node must be positive integer")
  }

} # end check_node_is_in_tree()

check_tree_is_valid <- function(tr){
  #`
  # Function description ------------------------------------------------------#
  # Test if a tree is valid- based on number of Nodes and Tips.
  #
  # Inputs:
  # tr: phylogenetic tree.
  #
  # Output:
  # None.
  #
  # Check input & function ----------------------------------------------------#
  if (class(tr) != 'phylo'){
    stop('Input must be a phylogenetic tree (object with class phylo)')
  }

  num_edges_for_node <- table(tr$edge)

  for (i in 1:Ntip(tr)){
    if (num_edges_for_node[i] != 1){
      stop(paste("Tip node", i, "has", num_edges_for_node[i], "edges. Should have 1 edge"))
    }
  }
  for (i in (Ntip(tr) + 1):(Nnode(tr) + Ntip(tr))){
    if (num_edges_for_node[i] != 2 && num_edges_for_node[i] != 3){
      stop(paste("Internal node", i, "has", num_edges_for_node[i], "edges. Should have 2(root) or 3 edge"))
    }
  }
}

check_convergence_possible <- function(vec, discrete_or_continuous){
  # Function description ------------------------------------------------------#
  # Test if vector, which represents values on the tips of tree, could plausibly
  # be consistent with convergence of those values. Eg. A value needs to appear
  # at leaset twice in the vector.
  #
  # Inputs:
  # vec. Vector. A vector of numbers. If "discrete" then vector must be binary.
  # discrete_or_continuous. Character. Either "discrete" or "continuous."
  #
  # Output:
  # None.
  #
  # Check input & function ----------------------------------------------------#
  if (discrete_or_continuous != "discrete"){
    if (discrete_or_continuous != "continuous"){
      stop("discrete_or_continuous must be a string 'discrete' or 'continuous.'")
    }
  }

  if (discrete_or_continuous == "discrete"){
    check_if_binary_vector(vec)
    if (sum(vec) >= (length(vec)-1) | sum(vec) <= 1){
      stop("Convergence is not possible for this phenotype")
    }

  } else if (discrete_or_continuous == "continuous"){
    if (length(unique(vec)) == 1){
      stop("Convergence is not possible for this phenotype")
    }
  }
} # end check_convergence_possible()

is_tip <- function(node_num, tr){
  # Function description ------------------------------------------------------#
  # Test if a node is a tip or an internal node.
  #
  # Inputs:
  # node_num: Integer. Index of node.
  # tr: phylogenetic tree.
  #
  # Output:
  # Logical. TRUE OR FALSE.
  #
  # Check input ---------------------------------------------------------------#
  check_tree_is_valid(tr)
  check_is_number(node_num)
  if (node_num < 1 ||node_num %% 1 != 0){
    stop("Node number must be a positive integer")
  }
  check_node_is_in_tree(node_num, tr)
  #
  # Function & return output --------------------------------------------------#
  return(node_num <= Ntip(tr))
} # end is_tip()


check_if_g_mat_can_be_plotted <- function(geno_matrix){
  # Function description -------------------------------------------------------
  # TODO
  # In order to make a heatmap there need to be 1) at least two columns, 2)
  # two different values within the matrix (0 and 1).
  # There can be NAs in the matrix.
  #
  # Inputs:
  # geno_matrix. Matrix. 1, 0, and/or NA.
  #
  # Outputs:
  # plot_logical. Logical. TRUE or FALSE.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(geno_matrix, min_rows = 1, min_cols = 2)

  if (sum(as.vector(geno_matrix)[!is.na(as.vector(geno_matrix))] %% 1 != 0) != 0){
    stop("Joint genotype matrix + phenotype must contain only 1, 0, or NA. (For discrete heatmap plot).")
  }

  # Function -------------------------------------------------------------------
  ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
  zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
  nas <- sum(is.na(geno_matrix)) > 1

  plot_logical <- FALSE #
  if (ones == 1 && zeros == 1 && nas == 0) {
    plot_logical <- TRUE
  }
  if (ones + zeros + nas == 3) {
    plot_logical <- TRUE
  }

  # Return output --------------------------------------------------------------
  if(!is.logical(plot_logical)){stop("Output must be a logical")}
  return(plot_logical)
} # end check_if_g_mat_can_be_plotted()


check_str_is_discrete_or_continuous <- function(input){
  # Function description -------------------------------------------------------
  # Check if the string is either "discrete" or "continuous"
  #
  # Inputs:
  # input. String. "discrete" or "continuous"
  #
  # Outputs:
  # None.
  #
  # Check input ----------------------------------------------------------------
  check_is_string(input)

  # Function -------------------------------------------------------------------
  if (input != "discrete"){
    if (input != "continuous"){
      stop("Input must be either 'discrete' or 'continuous'.")
    }
  }
} # end check_str_is_discrete_or_continuous()

check_if_ancestral_reconstruction_method_compatible_with_ape <- function(input){
  # Function description -------------------------------------------------------
  # Check that the reconstruction method that is being fed to ape::ace() is
  # one of the four acceptable methods.
  #
  # Inputs:
  # input. String.
  #
  # Output:
  # none.

  # Check inputs -------------------------------------------------------------
  check_is_string(input)

  # Function -----------------------------------------------------------------
  acceptable_methods <- c("ML", "REML", "pic", "GLS")
  if (!input %in% acceptable_methods){
    stop("Reconstruction methods for ape::ace must be either: ML, REML, pic, or GLS.")
  }
} # end check_if_ancestral_reconstruction_method_compatible_with_ape()
