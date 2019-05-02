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
  if (a >= 1 | a <= 0){
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

check_input_format <- function(pheno, tr, geno, name, dir, perm, a, annot){
  # Function description -------------------------------------------------------
  # Check that all of the inputs into phyC are in the valid format.
  #
  # Input:
  # pheno. Matrix. Phenotype.
  # tr.    Phylo. Tree.
  # geno.  Matrix. Genotype.
  # name.  Character. Output name.
  # dir.   Character. Output path.
  # perm.  Number. Times to shuffle the data on the tree to create a null distribution for the permutation test.
  # a.     Number. Alpha.
  # annot. Matrix. Annotation for heatmaps.

  # Output:
  # discrete_or_continuous. Character. Either "discrete" or "continuous". Describes the input phenotype.
  #
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, Ntip(tr), 2, NULL, 2) # Genotype matrix should have same rows as tr$tip.label, at least 2 genotypes in the columns
  check_dimensions(pheno, Ntip(tr), 2, 1, 1) # Phnoetype matrix should have same rows as tr$tip.label and exactly 1 column
  check_rownames(geno, tr) # Genotype rownames should be in the same order as the tr$tip.label
  check_rownames(pheno, tr) # Phenotype rownames should be in the same order as the tr$tip.label
  check_for_NA_and_inf(geno)
  check_for_NA_and_inf(pheno)
  check_for_root_and_boostrap(tr)
  check_tree_is_valid(tr)
  check_if_binary_matrix(geno)
  check_is_string(name)
  check_if_dir_exists(dir)
  check_if_permutation_num_valid(perm)
  check_if_alpha_valid(a)
  if (!is.null(annot)){
    check_dimensions(annot, Ntip(tr), 2, 2, 2)
  }

  # Function -------------------------------------------------------------------
  discrete_or_continuous <- assign_pheno_type(pheno)

  # Check and return output ----------------------------------------------------
  check_is_string(discrete_or_continuous)
  return(discrete_or_continuous)
} # end check_input_format()

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

check_for_root_and_boostrap <- function(tr){
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
} # end check_for_root_and_boostrap()

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
  if (sum(!(vec %in% c(0, 1))) > 0 | class(vec) != "integer"){
    stop("Vector should be only 1s and 0s")
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
  if (sum(row.names(mat) != tr$tip.label) != 0){
    stop("Matrix must be formatted with samples in matrix in the same order as tree$tip.label.")
  }
} # end check_rownames()

check_is_number <- function(num){
  # Function description -------------------------------------------------------
  # Check that input is some type of number.
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
  check_for_root_and_boostrap(tr)
  check_is_number(node_val)

  if (node_val > Nnode(tr) + Ntip(tr)){
    stop("Node number is too high; not found in tree.")
  }
  if (node_val < 1 | !is.integer(node_val)){
    stop("Node number must be positive integer")
  }
} # end check_node_is_in_tree()

check_tree_is_valid <- function(tr){
  #` A valid tree has N nodes, with n_tips nodes being tip nodes, numbered 1 through n_tips`

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
  return(TRUE)
}

# check_transitions <- function(transition_vector, tr){
#   for (i in 1:length(transition_vector)){
#     parent <- tr$edge[i, 1]
#     child <- tr$edge[i, 2]
#
#   }
#
# } # end check_transitions()

check_convergence_possible <- function(vec, discrete_or_continuous){
  convergence_not_possible <- FALSE
  if (discrete_or_continuous == "discrete"){
    check_if_binary_vector(vec)
    if (sum(vec) >= (length(vec)-1) | sum(vec) <= 1){
      convergence_not_possible <- TRUE
    }

  } else if (discrete_or_continuous == "continuous"){
    if (length(unique(vec)) == 1){
      convergence_not_possible <- TRUE
    }
  }
  if (convergence_not_possible){
    stop("Convergence is not possible for this phenotype")
  }
} # end check_convergence_possible()
