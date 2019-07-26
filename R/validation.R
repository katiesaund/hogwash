#' check_dimensions
#'
#' @description Check that the main input is a matrix and that said matrix is of
#'  the specified dimensions.
#' @param mat Matrix.
#' @param exact_rows Numeric. Describes expected number of rows in matrix.
#'  Default is NULL.
#' @param min_rows Numeric. Describes minimum number of rows in matrix. Must be
#'  specified.
#' @param exact_cols Numeric. Describes expected number of columns in matrix.
#'  Default NULL.
#' @param min_cols Numeric. Describes minimum number of columns in matrix. Must
#'  be specified.
#'
#' @noRd
check_dimensions <- function(mat,
                             exact_rows = NULL,
                             min_rows,
                             exact_cols = NULL,
                             min_cols){
  # Check input ----------------------------------------------------------------
  if (!is.null(exact_rows)) {
    check_is_number(exact_rows)
  }
  check_is_number(min_rows)
  if (!is.null(exact_cols)) {
    check_is_number(exact_cols)
  }
  check_is_number(min_cols)
  if (class(mat) != "matrix") {
    stop("input must be a matrix")
  }

  # Function -------------------------------------------------------------------
  if (nrow(mat) < min_rows) {
    stop("matrix has too few rows")
  }
  if (ncol(mat) < min_cols) {
    stop("matrix has too few columns")
  }
  if (!(is.null(exact_rows))) {
    if (nrow(mat) != exact_rows) {
      stop("matrix has wrong number of rows")
    }
  }
  if (!(is.null(exact_cols))) {
    if (ncol(mat) != exact_cols) {
      stop("matrix has wrong number of columns")
    }
  }
} # end check_dimensions()

#' check_num_between_0_and_1
#'
#' @description Check that the input is within a the expected range
#'  (0 < num < 1). To be used to check values such as fdr, p-value, and
#'  bootstrap confidence.
#'
#' @param num Number.
#'
#' @noRd
check_num_between_0_and_1 <- function(num){
  # Check input ----------------------------------------------------------------
  check_is_number(num)

  # Function -------------------------------------------------------------------
  if (num > 1 | num < 0) {
    stop("Provide a value between 0 and 1.")
  }
} # end check_num_between_0_and_1()

#' check_if_dir_exists
#'
#' @description Check that output directory exists so data can be saved in it.
#'
#' @param dir Character. Path to output directory.
#'
#' @noRd
check_if_dir_exists <- function(dir){
  # Check input ----------------------------------------------------------------
  check_is_string(dir)

  # Function -------------------------------------------------------------------
  if (!dir.exists(dir)) {
    stop("The output directory indicated does not exist.")
  }
} # end check_if_dir_exists()

#' check_if_permutation_num_valid
#'
#' @description Check that the permutation number indicated is valid
#'  (1 <= perm). A typical choice is 10,000.
#'
#' @param perm Integer. Number of times to shuffle the data on the tree to
#'  create a null distribution for the permutation test.
#'
#' @noRd
check_if_permutation_num_valid <- function(perm){
  # Check input ----------------------------------------------------------------
  check_is_number(perm)

  # Function -------------------------------------------------------------------
  if (perm < 1 || perm %% 1 != 0) {
    stop("The permutation number should be a positive integer indicating the
         number of null distributions to create.")
  }
} # end check_if_permutation_num_valid()

#' check_is_string
#'
#' @description Check that the input is a character string.
#'
#' @param char Character.
#'
#' @noRd
check_is_string <- function(char){
  # Check input & function -----------------------------------------------------
  if (!is.character(char)) {
    stop("Object must be a character string.")
  }
} # end check_is_string()

#' check_if_vector
#'
#' @description Check that input is a vector.
#'
#' @param vector Vector.
#'
#' @noRd
check_if_vector <- function(vector){
  # Check input & function -----------------------------------------------------
  if (!is.vector(vector)) {
    stop("Input must be a vector")
  }
} # end check_if_vector()

#' check_for_NA_and_inf
#'
#' @description Check that matrix contains no NAs and no +/- infinities.
#'
#' @param mat Matrix.
#'
#' @noRd
check_for_NA_and_inf <- function(mat){
  # Check input & function -----------------------------------------------------
  if (class(mat) != "matrix") {
    stop("Input should be a matrix.")
  }
  if (sum(is.na(mat)) > 0) {
    stop("Input matrices should not have any NA values.")
  }
  if (sum(mat == -Inf) > 0) {
    stop("Inpute matrices should not have any -Inf values.")
  }
  if (sum(mat == Inf) > 0) {
    stop("Inpute matrices should not have any -Inf values.")
  }
} # end check_for_NA_and_inf()

#' check_for_root_and_bootstrap
#'
#' @description Check that phylogenetic tree is rooted and contains bootstrap
#'  values in the node labels. Trees need to have roots for the ancestral
#'  reconstruction function ape::ace() to work properly. Trees need to have
#'  bootstrap values so that confidence in the tree edges can be measured.
#'
#' @param tr Phylo.
#'
#' @noRd
check_for_root_and_bootstrap <- function(tr){
  # Check input & function -----------------------------------------------------
  if (class(tr) != "phylo") {
    stop("Tree must be phylo object")
  }
  if (!ape::is.rooted(tr)) {
    stop("Tree must be rooted")
  }
  if (is.null(tr$node.label)) {
    stop("Tree must have bootstrap values in the nodes")
  }
  if (length(tr$node.label) != ape::Nnode(tr)) {
    stop("Tree must have bootstrap values for each node")

  }
  for (i in 2:length(tr$node.label)) {
    check_is_number(tr$node.label[i]) # Node labels must be bootstrap values
    if (tr$node.label[i] < 0) {
      stop("Tree bootstrap values must be >= 0")
    }
  }
} # end check_for_root_and_bootstrap()

#' check_if_binary_vector
#'
#' @description Check that the matrix only contains values 1 or 0. Useful for
#'  checking things for binary phenotypes, confidence, $transition, etc...
#'
#' @param vec Vector
#'
#' @noRd
check_if_binary_vector <- function(vec) {
  # Check input & function -----------------------------------------------------
  if (sum(!(vec %in% c(0, 1))) > 0) {
    stop("Vector should be only 1s and 0s")
  }
  if (class(vec) != "integer") {
    if (class(vec) != "numeric") {
      stop("Vector should be only 1s and 0s")
    }
  }
} # end check_if_binary_vector()

#' check_if_binary_vector_numeric
#'
#' @description Check that the matrix only contains values 1 or 0.
#'
#' @param vec Vector.
#'
#' @noRd
check_if_binary_vector_numeric <- function(vec) {
  # Check input & function -----------------------------------------------------
  if (sum(!(vec %in% c(0, 1))) > 0 | class(vec) != "numeric") {
    stop("Vector should be only 1s and 0s")
  }
} # end check_if_binary_vector_numeric()

#' check_if_binary_matrix
#'
#' @description Check that the matrix only contains values 1 or 0. Useful for
#'  binary phenotype matrix or genotype matrix.
#'
#' @param mat Matrix.
#'
#' @noRd
check_if_binary_matrix <- function(mat) {
  # Check input & function -----------------------------------------------------
  if (is.null(dim(mat))) {
    stop("Must be a matrix")
  }
  if (sum(!(mat %in% c(0, 1))) > 0 | class(mat) != "matrix") {
    stop("Genotype matrix should be only 1s and 0s")
  }
  if (nrow(mat) == 0 | ncol(mat) == 0) {
    stop("matrix must have columns and rows.")
  }
} # end check_if_binary_matrix()

#' check_file_exists
#'
#' @description Check that the file exists.
#'
#' @param file_name Character.
#'
#' @noRd
check_file_exists <- function(file_name) {
  # Check input & function -----------------------------------------------------
  if (!file.exists(file_name)) {
    stop("File does not exist")
  }
} # end check_file_exists()

#' check_rownames
#'
#' @description Check that phylogenetic tree tip labels are identical to the
#'  matrix row.names.
#'
#' @param mat Matrix. Nrow(mat) == Ntip(tr).
#' @param tr Phylo.
#'
#' @noRd
check_rownames <- function(mat, tr){
  # Check input ----------------------------------------------------------------
  if (class(mat) != "matrix" | class(tr) != "phylo") {
    stop("Inputs are incorrectly formatted.")
  }
  check_dimensions(mat,
                   exact_rows = ape::Ntip(tr),
                   min_rows = ape::Ntip(tr),
                   exact_cols = NULL,
                   min_cols = 1)

  # Function -------------------------------------------------------------------
  if (is.null(row.names(mat))) {
    stop("Matrix must have row.names()")
  }
  if (is.null(tr$tip.label))  {
    stop("Tree must have tip labels")
  }
  if (sum(row.names(mat) != tr$tip.label) != 0) {
    stop("Matrix must be formatted with samples in matrix in the same order as
         tree$tip.label.")
  }
} # end check_rownames()

#' check_is_number
#'
#' @description Check that input is some type of number. Check that the input is
#'  a single number.
#'
#' @param num Number. Could be numeric, double, or integer.
#'
#' @noRd
check_is_number <- function(num){
  # Check input & function -----------------------------------------------------
  if (!is.numeric(num)) {
    if (!is.integer(num)) {
      if (!is.double(num)) {
        stop("Must be a number")
      }
    }
  }
  if (is.data.frame(num)) {
    stop("Number can't be a dataframe")
  }
  if (is.matrix(num)) {
    stop("Number can't be a matrix")
  }
  if (length(num) != 1) {
    stop("Must be a single number")
  }
  if (num == Inf | num == -Inf) {
    stop("This can't handle infinity.")
  }
} # end check_is_number()

#' check_node_is_in_tree
#'
#' @description Test if a node value is plausibly contained within the tree.
#'
#' @param node_val Integer. Index of node.
#' @param tr Phylo.
#'
#' @noRd
check_node_is_in_tree <- function(node_val, tr){
  # Check input & function -----------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_is_number(node_val)

  if (node_val > ape::Nnode(tr) + ape::Ntip(tr)) {
    stop("Node number is too high; not found in tree.")
  }
  if (node_val < 1 | node_val %% 1 != 0) {
    stop("Node must be positive integer")
  }
} # end check_node_is_in_tree()

#' check_tree_is_valid
#'
#' @description Test if a tree has valid structure. Each tree node should touch
#'  either three edges (internal node) or two edges (root node). This function
#'  checks if a phylogenetic tree is structured correctly so that later
#'  functions will work as expected.
#'
#' @param tr Phylo.
#'
#' @noRd
check_tree_is_valid <- function(tr){
  # Check input & function -----------------------------------------------------
  if (class(tr) != "phylo") {
    stop("Input must be a phylogenetic tree (object with class phylo)")
  }

  num_edges_for_node <- table(tr$edge)

  for (i in 1:ape::Ntip(tr)) {
    if (num_edges_for_node[i] != 1) {
      stop(paste("Tip node", i, "has", num_edges_for_node[i],
                 "edges. Should have 1 edge"))
    }
  }
  for (i in (ape::Ntip(tr) + 1):(ape::Nnode(tr) + ape::Ntip(tr))) {
    if (num_edges_for_node[i] != 2 && num_edges_for_node[i] != 3) {
      stop(paste("Internal node", i, "has", num_edges_for_node[i],
                 "edges. Should have 2(root) or 3 edge"))
    }
  }
}

#' check_convergence_possible
#'
#' @description Test if vector, which represents values on the tips of tree,
#'  could plausibly be consistent with convergence of those values, e.g. A value
#'  needs to appear at leaset twice in the vector.
#'
#' @param vec Vector. A vector of numbers. If "discrete" then vector must be
#'  binary.
#' @param discrete_or_continuous Character. Either "discrete" or "continuous."
#'
#' @noRd
check_convergence_possible <- function(vec, discrete_or_continuous){
  # Check input & function -----------------------------------------------------
  if (discrete_or_continuous != "discrete") {
    if (discrete_or_continuous != "continuous") {
      stop("discrete_or_continuous must be a string 'discrete' or
           'continuous.'")
    }
  }

  if (discrete_or_continuous == "discrete") {
    check_if_binary_vector(vec)
    if (sum(vec) >= (length(vec) - 1) | sum(vec) <= 1) {
      stop("Convergence is not possible for this phenotype")
    }

  } else if (discrete_or_continuous == "continuous") {
    if (length(unique(vec)) == 1) {
      stop("Convergence is not possible for this phenotype")
    }
  }
} # end check_convergence_possible()

#' is_tip
#'
#' @description Test if a node is a tree tip. An internal node should return
#'  false.
#'
#' @param node_num Integer. Index of node.
#' @param tr Phylo.
#'
#' @return Logical. TRUE OR FALSE.
#' @noRd
is_tip <- function(node_num, tr){
  # Check input ----------------------------------------------------------------
  check_tree_is_valid(tr)
  check_is_number(node_num)
  if (node_num < 1 || node_num %% 1 != 0) {
    stop("Node number must be a positive integer")
  }
  check_node_is_in_tree(node_num, tr)

  # Function & return output ---------------------------------------------------
  return(node_num <= ape::Ntip(tr))
} # end is_tip()


#' check_if_g_mat_can_be_plotted
#'
#' @description The program cannot plot all results. In order to plot the
#'  heatmap results there needs to be 1) at least two columns in the matrix, 2)
#'  two different values within the matrix (0 and 1). There can be NAs in the
#'  matrix.
#'
#' @param geno_matrix Matrix. 1, 0, and/or NA.
#'
#' @return plot_logical. Logical. TRUE or FALSE.
#' @noRd
check_if_g_mat_can_be_plotted <- function(geno_matrix){
  # Check input & function -----------------------------------------------------
  if (sum(class(geno_matrix) != "data.frame",
          class(geno_matrix) != "matrix") == 2) {
    # Neither matrix nor dataframe
    plot_logical <- FALSE
  } else {
    # Either a matrix or dataframe
    if (nrow(geno_matrix) < 1 | ncol(geno_matrix) < 2) {
      # Matrix/dataframe is too small for heatmap to plot
      plot_logical <- FALSE
    } else {
      # Matrix/dataframe is big enough for heatmap to plot
      if (sum(as.vector(geno_matrix)[!is.na(as.vector(geno_matrix))] %% 1 != 0)
          != 0) {
        # Matrix/dataframe contains invalid values
        stop("Joint genotype matrix + phenotype must contain only 1, 0, or NA.
             (For discrete heatmap plot).")
      }
      # Matrix/dataframe contains only valid values
      ones <- sum(geno_matrix == 1, na.rm = TRUE) > 1
      zeros <- sum(geno_matrix == 0, na.rm = TRUE) > 1
      nas <- sum(is.na(geno_matrix)) > 1
      plot_logical <- FALSE
      if (ones == 1 && zeros == 1 && nas == 0) {
        # Has at least two values to plot
        plot_logical <- TRUE
      } else if (ones + zeros + nas == 3) {
        # Has at least two values to plot
        plot_logical <- TRUE
      } else {
        # Does not have at least two values to plot
        plot_logical <- FALSE
      }
    }
  }
  # Return output --------------------------------------------------------------
  return(plot_logical)
} # end check_if_g_mat_can_be_plotted()

#' check_str_is_discrete_or_continuous
#'
#' @description Check if the string is either "discrete" or "continuous."
#'
#' @param input String. "discrete" or "continuous"
#'
#' @noRd
check_str_is_discrete_or_continuous <- function(input){
  # Check input ----------------------------------------------------------------
  check_is_string(input)

  # Function -------------------------------------------------------------------
  if (input != "discrete") {
    if (input != "continuous") {
      stop("Input must be either 'discrete' or 'continuous'.")
    }
  }
} # end check_str_is_discrete_or_continuous()

#' check_equal
#'
#' @description This function simply replaces the if(){stop()} statements for
#'  asserting that two things are equal. Example:
#'   check_equal(nrow(genotype_matrix), Ntip(tree)) should not give any
#'   warnings.
#'
#' @param first_number Number.
#' @param second_number Number.
#'
#' @noRd
check_equal <- function(first_number, second_number){
  # Check inputs ---------------------------------------------------------------
  check_is_number(first_number)
  check_is_number(second_number)

  # Function -------------------------------------------------------------------
  if (first_number != second_number) {
    stop("Inputs are not equal")
  }
} # end check_equal()

#' check_anc_rec_compatible
#'
#' @description Check that the reconstruction method that is being fed to
#'  ape::ace() is one of the four acceptable methods. The four methods are:
#'  "ML", "REML", "pic", and "GLS." For the intial implementation of this
#'  package the default (hard-coded) option is always maximum likelihood ("ML").
#' @param input String. Either "ML", "REML", "pic", or "GLS."
#'
#' @noRd
check_anc_rec_compatible <- function(input){
  # Check inputs -------------------------------------------------------------
  check_is_string(input)
  check_equal(length(input), 1)

  # Function -----------------------------------------------------------------
  acceptable_methods <- c("ML", "REML", "pic", "GLS")
  if (!input %in% acceptable_methods) {
    stop("Reconstruction methods for ape::ace must be either:
         ML, REML, pic, or GLS.")
  }
} # end check_anc_rec_compatible()

#' check_class
#'
#' @param obj  Any R object.
#' @param cls Character string. Describes a class type, e.g. "matrix", "list",
#'  "vector", etc...
#'
#' @noRd
check_class <- function(obj, cls){
  # Check inputs -------------------------------------------------------------
  check_is_string(cls)

  # Function -----------------------------------------------------------------
  if (class(obj) != cls) {
    stop("Object does not have expected class.")
  }
} # end check_class()

#' check_str_is_test_name
#'
#' @description Verify that the string is one of the three valid gwas test
#'  names:
#'  * "continuous"
#'  * "synchronous"
#'  * "phyc"
#'  If not valid, give error message.
#'
#' @param test_name Character.
#'
#' @noRd
check_str_is_test_name <- function(test_name){
  # Check input ----------------------------------------------------------------
  check_class(test_name, "character")

  # Function -------------------------------------------------------------------
  valid_names <- c("continuous", "synchronous", "phyc")
  if (!(test_name %in% valid_names)) {
    stop("Test name must be either continuous, synchronous, or phyc")
  }
} # end check_str_is_test_name()
