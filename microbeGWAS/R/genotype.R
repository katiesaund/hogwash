
get_dropped_genotypes <- function(geno, keepers){
  # Function description -------------------------------------------------------
  # TODO
  # Compute ancestral reconstruction from a continuous or discrete input.
  #
  # Inputs:
  # Varname. Var class. Description.
  #
  # Outputs:
  # "anc_rec" = ML_anc_rec. Vector. Description.
  #
  # Check input ----------------------------------------------------------------

  # Function -------------------------------------------------------------------

  # Return output --------------------------------------------------------------
  dropped_genotype_names <- colnames(geno)[!keepers]
  return(dropped_genotype_names)
} # end get_dropped_genotypes



reduce_redundancy <- function(mat, tr, dir, name){
  # Function description -------------------------------------------------------
  # This function REMOVES GENOTYPES THAT ARE: RARE, TOO COMMON, OR IDENTICAL
  #
  # Inputs:
  # mat.  Character. Path to matrix file.
  # tr.   Phylo.
  # dir.  Character.
  # name. Character.
  #
  # Output:
  # mat. Matrix.
  #
  # Check inputs ---------------------------------------------------------------
  check_for_root_and_bootstrap(tr)
  check_if_dir_exists(dir)
  check_is_string(name)
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
