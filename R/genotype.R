#' Get names of genotypes dropped from analysis
#'
#' @description  Given a matrix with both desirable genotypes (keepers) and
#'  undesirable genotypes and a numeric index of the keepers, get the names of
#'  the undesirable genotypes.
#'
#' @param geno Matrix. Numeric, binary genotype matrix. Columns = genotypes.
#'   Rows = samples.
#' @param keepers Logical vector. Length == length(genotype_transition)  ==
#'  length(genotype_confidence) == Number of genotypes. True indicates the
#'  genotype has at least 2 high confidence genotype transition edges. False
#'  indicates genotype lacks at least 2 high confidence genotype transition
#'  edges.
#'
#' @return dropped_genotype_names. Character vector. Names of the non-keepers
#'  (genotypes not to be processed downstream).
#'
#' @noRd
#'
get_dropped_genotypes <- function(geno, keepers){
  # Check input ----------------------------------------------------------------
  check_if_binary_matrix(geno)
  check_class(keepers, "logical")
  check_if_vector(keepers)
  if (ncol(geno) != length(keepers)) {
    stop("Keepers must have an entry for each genotype.")
  }
  if (sum(keepers) == 0) {
    stop("There are no genotypes left to test.")
  }

  # Function -------------------------------------------------------------------
  dropped_genotype_names <- colnames(geno)[!keepers]

  # Return output --------------------------------------------------------------
  return(dropped_genotype_names)
} # end get_dropped_genotypes()

#' Remove rare or common genotypes from genotype matrix
#'
#' @description This function removes genotypes are either too rare (only occur
#'  once) or too common (occur everywhere or everywhere but once). This step is
#'  necessary so that ape::ace() doesn't break and so that only genotypes where
#'  convergence is possible are included for testing.
#'
#' @param mat Matrix. Columns correspond to genotypes. Rows correspond to tree
#'  tips. Binary.
#' @param tr Phylo. Tree.
#'
#' @return A list with the following items:
#'   \describe{
#'     \item{mat}{Matrix.}
#'     \item{dropped_genotype_names}{Character vector.}
#'   }
#' @noRd
#'
remove_rare_or_common_geno <- function(mat, tr){
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(mat)

  # Function -------------------------------------------------------------------
  geno_to_drop <- rep(FALSE, ncol(mat))
  geno_to_drop[colSums(mat) <= 1] <- TRUE
  geno_to_drop[colSums(mat) == (ape::Ntip(tr) - 1)] <- TRUE
  geno_to_drop[colSums(mat) == ape::Ntip(tr)] <- TRUE
  dropped_genotype_names <- colnames(mat)[geno_to_drop]
  mat <- mat[, !geno_to_drop, drop = FALSE]

  # Check and return output ----------------------------------------------------
  check_if_binary_matrix(mat)
  check_rownames(mat, tr)
  results <- list("mat" = mat,
                  "dropped_genotype_names" = dropped_genotype_names)
  return(results)
} # end remove_rare_or_common_geno()
