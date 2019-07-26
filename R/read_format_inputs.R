#' check_input_format
#'
#' @description Check that all of the inputs into hogwash are in the valid
#'  format.
#' @param pheno  Matrix. Phenotype. Phenotype. matrix should have same rows as
#'  tr$tip.label and exactly 1 column. Phenotype rownames should be in the same
#'  order as the tr$tip.label.
#' @param tr Phylo. Tree.
#' @param geno Matrix. Genotype. Genotype matrix should have same rows as
#' tr$tip.label, at least 2 genotypes in the columns. Genotype rownames should
#'  be in the same order as the tr$tip.label.
#' @param name Character. Output name.
#' @param dir  Character. Output path.
#' @param perm  Number. Times to shuffle the data on the tree to create a null
#'  distribution for the permutation test.
#' @param fdr Number. False discovery rate. Between 0 and 1.
#'
#' @return discrete_or_continuous. Character. Either "discrete" or "continuous".
#'  Describes the input phenotype.
#' @noRd
check_input_format <- function(pheno, tr, geno, name, dir, perm, fdr){
  # Check input ----------------------------------------------------------------
  check_dimensions(geno, ape::Ntip(tr), 2, NULL, 2)
  check_dimensions(pheno, ape::Ntip(tr), 2, 1, 1)
  check_rownames(geno, tr)
  check_rownames(pheno, tr)
  check_for_NA_and_inf(geno)
  check_for_NA_and_inf(pheno)
  check_if_binary_matrix(geno)
  check_is_string(name)
  check_if_dir_exists(dir)
  check_if_permutation_num_valid(perm)
  check_num_between_0_and_1(fdr)

  # Function -------------------------------------------------------------------
  discrete_or_continuous <- assign_pheno_type(pheno)

  # Check and return output ----------------------------------------------------
  check_is_string(discrete_or_continuous)
  return(discrete_or_continuous)
} # end check_input_format()
