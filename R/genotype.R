get_dropped_genotypes <- function(geno, keepers){
  # Function description -------------------------------------------------------
  # Given a matrix with both desirable genotypes (keepers) and undesirable
  # genotypes and a numeric index of the keepers, get the names of the keepers.
  #
  # Inputs:
  # geno. Matrix. Numeric, binary matrix.
  # keepers. ?
  #
  # Outputs:
  # dropped_genotype_names. Character vector. Names of the non-keepers
  #                         (genotypes not to be processed downstream).
  #
  # Check input ----------------------------------------------------------------
  check_if_binary_matrix(geno)
  # TODO add stuff for keepers

  # Function -------------------------------------------------------------------
  dropped_genotype_names <- colnames(geno)[!keepers]

  # Return output --------------------------------------------------------------
  return(dropped_genotype_names)
} # end get_dropped_genotypes

reduce_redundancy <- function(mat, tr){
  # Function description -------------------------------------------------------
  # This function removes genotypes are either too rare (only occur once) or too
  # common (occur everywhere or everywhere but once).
  #
  # Inputs:
  # mat.  Character. Path to matrix file.
  # tr.   Phylo.
  #
  # Output:
  # mat. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_binary_matrix(mat)

  # Function -------------------------------------------------------------------
  geno_to_drop <- rep(FALSE, ncol(mat))
  geno_to_drop[colSums(mat) <= 1] <- TRUE
  geno_to_drop[colSums(mat) == (Ntip(tr) - 1)] <- TRUE
  geno_to_drop[colSums(mat) == Ntip(tr)] <- TRUE
  dropped_genotype_names <- colnames(mat)[geno_to_drop]
  mat <- mat[ , !geno_to_drop, drop = FALSE]

  # Check and return output ----------------------------------------------------
  check_if_binary_matrix(mat)
  check_rownames(mat, tr)
  results <- list("mat" = mat, "dropped_genotype_names" = dropped_genotype_names)
  return(results)
} # end reduce_redundancy()
