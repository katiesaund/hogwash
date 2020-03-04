#' Format phylogenetic tree
#'
#' @description For a phylogenetic tree convert all bootstrap support values to
#'   numeric (from character). Typically, a tree's root doesn't have a value but
#'   this adds a zero to that location (if necessary) because functions later in
#'   hogwash require that all nodes have numeric values.
#'
#' @param tr Phylo.
#'
#' @return tr. Phylo.
#'
#' @noRd
format_tree <- function(tr){
  # Check input ----------------------------------------------------------------
  check_class(tr, "phylo")
  if (is.null(tr$node.label)) {
    stop("trees must have support values at the nodes")
  }
  if (!ape::is.rooted(tr)) {
    tr <- phytools::midpoint.root(tr)
  }

  # Function -------------------------------------------------------------------
  for (i in 1:ape::Nnode(tr)) {
    if (tr$node.label[i] == "") {
      tr$node.label[i] <- 0
    } else if (tr$node.label[i] == "Root") {
      tr$node.label[i] <- 0
    }
  }
  tr$node.label <- as.numeric(tr$node.label)

  # Check and return output ----------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_tree_is_valid(tr)
  return(tr)
}

#' Convert matrix to a vector
#'
#' @description Convert a single column matrix into a vector, retain row names
#'   as names of vector.
#'
#' @param mat Matrix. Matrix should have only 1 column.
#'
#' @return vec. Named vector.
#'
#' @noRd
convert_matrix_to_vector <- function(mat){
  # Check input ----------------------------------------------------------------
  check_dimensions(mat, NULL, 1, 1, 1)

  # Function -------------------------------------------------------------------
  vec <- as.vector(unlist(mat[, 1]))
  names(vec) <- row.names(mat)

  # Check and return output ----------------------------------------------------
  check_if_vector(vec)
  return(vec)
}

#' Convert phenotype matrix to vector
#'
#' @description Prepare phenotype for downstream analyses. Convert phenotype
#'   matrix to phenotype vector. Check if the continuous phenotype is normally
#'   distributed.
#'
#' @param pheno Matrix with 1 column.
#' @param disc_cont String. Either "discrete" or "continuous."
#' @param tr Phylo. Ntip = nrow(pheno).
#'
#' @return pheno_vector. Vector. Length = Ntip(tr).
#'
#' @noRd
prepare_phenotype <- function(pheno, disc_cont, tr){
  # Check input ----------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_dimensions(pheno,
                   exact_rows = ape::Ntip(tr),
                   exact_cols = 1,
                   min_rows = 1,
                   min_cols = 1)
  check_str_is_discrete_or_continuous(disc_cont)

  # Function -------------------------------------------------------------------
  check_if_phenotype_normal(pheno, disc_cont)
  pheno_vector <- convert_matrix_to_vector(pheno)
  check_convergence_possible(pheno_vector, disc_cont)
  internal_report_phylogenetic_signal(pheno, tr)

  # Check and return output ----------------------------------------------------
  check_equal(length(pheno_vector), ape::Ntip(tr))
  return(pheno_vector)
}
